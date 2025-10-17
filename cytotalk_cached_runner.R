## ============================================================
## CytoTalk cached runner with FILE-BACKED MI (HDF5)
##  - Avoids RAM blowups by writing MI to disk (per cell type)
##  - Two ARACNE modes:
##      * "full": load MI fully (exactly matches parmigene::aracne.m)
##      * "stream": blockwise DPI pruning without loading all of MI at once
##        (default; set topk_per_node = Inf for full-neighborhood DPI)
## ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(doParallel)
  library(foreach)
  library(parmigene)
  library(CytoTalk)
  library(rhdf5)
})

## --------- Progress tick ----------
.ct_tick <- local({
  t0 <- Sys.time(); i <- 0L
  function(step_msg) {
    i <<- i + 1L
    dt <- format(Sys.time() - t0, digits = 4)
    cat(sprintf("[%d / 8] (%s) %s\n", i, dt, step_msg)); flush.console()
  }
})

.ct_zero_diag <- function(m) { diag(m) <- 0; m }

## --------- Discretization (ensure numeric) ----------
.ct_discretize_dense <- function(mat_gxc) {
  mat_disc <- CytoTalk::discretize_sparse(Matrix::t(mat_gxc))  # cells x genes (sparse)
  if (storage.mode(mat_disc@x) != "double") storage.mode(mat_disc@x) <- "double"
  dense <- as.matrix(mat_disc)                                  # dense double
  storage.mode(dense) <- "double"
  if (is.null(colnames(dense))) colnames(dense) <- rownames(mat_gxc)
  dense
}

## --------- Expression filter (replace unexported subset_non_zero_old) ----------
.expr_filter <- function(mat, cutoff) {
  ncells <- ncol(mat)
  keep <- Matrix::rowSums(mat > 0) >= (cutoff * ncells)
  mat[keep, , drop = FALSE]
}

## ============================================================
## HDF5-backed MI utilities
## ============================================================

# Create an HDF5 dataset for MI (square p x p, doubles)
.h5_create_mi <- function(h5_path, dset, p) {
  if (file.exists(h5_path)) file.remove(h5_path)
  h5createFile(h5_path)
  h5createDataset(file = h5_path, dataset = dset, dims = c(p, p), storage.mode = "double", chunk = c(min(p, 512), min(p, 512)), level = 6)
}

# Write a symmetric block to HDF5
.h5_write_sym_block <- function(h5_path, dset, iblk, jblk, block_mat) {
  # block_mat: |iblk| x |jblk|
  h5write(block_mat, h5_path, dset, index = list(iblk, jblk))
  if (!identical(iblk, jblk)) {
    h5write(t(block_mat), h5_path, dset, index = list(jblk, iblk))
  }
}

# Read a row (as numeric vector) from HDF5
.h5_read_row <- function(h5_path, dset, i, p) {
  as.numeric(h5read(h5_path, dset, index = list(i, 1:p)))
}

# Read a rectangular block
.h5_read_block <- function(h5_path, dset, iblk, jblk) {
  h5read(h5_path, dset, index = list(iblk, jblk))
}

# Compute and store MI in HDF5 in blocks
.compute_mi_to_hdf5 <- function(dense_disc, h5_path, dset = "MI", block_size = 512, method_mi = "mm", ncores = NULL) {
  p <- ncol(dense_disc)  # genes
  .h5_create_mi(h5_path, dset, p)
  
  # Precompute entropy for all columns
  # infotheo::entropy expects a vector; dense_disc is cells x genes
  ent <- apply(dense_disc, 2, infotheo::entropy, method = method_mi)
  
  # Parallel backend controls foreach inside our own loops if we use it
  if (is.null(ncores) || ncores < 1L) ncores <- max(1L, parallel::detectCores() - 1L)
  doParallel::registerDoParallel(ncores)
  on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  
  # Blockwise computation of upper triangle (including diag)
  breaks <- seq.int(1L, p, by = block_size)
  for (bi in breaks) {
    ii <- bi:min(bi + block_size - 1L, p)
    # pre-extract X block for speed
    X <- dense_disc[, ii, drop = FALSE]
    h1 <- ent[ii]
    for (bj in breaks) {
      if (bj < bi) next  # upper triangle only
      jj <- bj:min(bj + block_size - 1L, p)
      Y <- dense_disc[, jj, drop = FALSE]
      h2 <- ent[jj]
      # compute MI block: for each i in ii, j in jj
      # I(X;Y) = H(X) + H(Y) - H(X,Y)
      # We'll compute H(X,Y) row-wise using infotheo::entropy on 2-col data.frame
      block <- matrix(0, nrow = length(ii), ncol = length(jj))
      for (a in seq_along(ii)) {
        x <- as.numeric(X[, a])
        ha <- h1[a]
        # vectorized over jj won't work directly; loop jj
        block[a, ] <- vapply(seq_along(jj), function(b) {
          y <- as.numeric(Y[, b])
          hb <- h2[b]
          hm <- min(ha, hb)
          hxy <- infotheo::entropy(data.frame(x, y), method = method_mi)
          if (hm == 0) 0 else ( (ha + hb) - hxy )
        }, numeric(1))
      }
      # zero diag if square block includes diagonal
      if (bi == bj) {
        diag(block) <- 0
      }
      .h5_write_sym_block(h5_path, dset, ii, jj, block)
      rm(X, Y, block); gc(FALSE)
    }
  }
  invisible(p)
}

## ============================================================
## ARACNE (DPI) with HDF5 MI
## ============================================================
# Two modes:
#  - "full": read MI fully into RAM and call parmigene::aracne.m (identical)
#  - "stream": perform DPI pruning by streaming MI rows; optional topk_per_node
#              (set topk_per_node = Inf to consider all neighbors for DPI)
#
# Returns a sparse adjacency matrix (dgCMatrix), non-negative, symmetric.
#
.aracne_from_hdf5 <- function(h5_path, dset = "MI", mode = c("stream", "full"), topk_per_node = Inf, dpi_tolerance = 0) {
  mode <- match.arg(mode)
  # read dims by reading small slice
  # Better: h5readAttributes to get dims (rhdf5 doesn't store dims attr here); we pass p down
  # We'll infer p as the maximum dimension length from the dataset listing:
  info <- h5dump(h5_path, load = FALSE)
  if (!dset %in% names(info)) stop("Dataset '", dset, "' not found in HDF5 file.")
  # Unfortunately h5dump doesn't expose dims without load; we'll read 1x1 and use dim():
  one <- h5read(h5_path, dset, index = list(1:1, 1:1))
  # Now read a row to get ncol:
  test_row <- h5read(h5_path, dset, index = list(1:1, NULL))
  p <- length(test_row); rm(one, test_row)
  
  if (mode == "full") {
    # WARNING: will load full MI into RAM (size ~ p^2 * 8 bytes)
    MI <- h5read(h5_path, dset)
    MI <- .ct_zero_diag(MI)
    adj <- Matrix::Matrix(parmigene::aracne.m(MI), sparse = TRUE)
    return(adj)
  }
  
  # STREAM MODE:
  # Build adjacency by DPI pruning using streaming rows.
  # Outline:
  # 1) For each node i, get MI(i, :)
  # 2) Choose neighbors N(i): all j with MI>0 (or topk_per_node largest)
  # 3) Initialize edges (i,j) for j in N(i)
  # 4) Apply DPI on (i,j,k) triangles by streaming: for each pair (j,k) in N(i),
  #    if MI(j,k) >= min(MI(i,j), MI(i,k)) - tol, drop the weaker of (i,j) or (i,k).
  #    We fetch MI(j,k) on demand from HDF5 (no big RAM)
  #
  # Complexity remains high but memory stays bounded. For speed, you can set a finite topk_per_node.
  
  # We'll build i<j edges only and symmetrize at the end
  keep_i <- vector("list", p)
  
  for (i in seq_len(p)) {
    mi_i <- .h5_read_row(h5_path, dset, i, p)  # length p
    mi_i[i] <- 0
    # neighbor candidates
    idx_nz <- which(mi_i > 0)
    if (length(idx_nz) == 0L) { keep_i[[i]] <- integer(0); next }
    # top-k if requested
    if (is.finite(topk_per_node) && length(idx_nz) > topk_per_node) {
      ord <- order(mi_i[idx_nz], decreasing = TRUE)
      idx_nz <- idx_nz[ord[seq_len(topk_per_node)]]
    }
    # start with all neighbors as kept
    kept <- rep(TRUE, length(idx_nz))
    
    # DPI pruning against triangles with center i
    # For all pairs (j,k) among neighbors:
    if (length(idx_nz) >= 2L) {
      for (a in seq_along(idx_nz)) {
        j <- idx_nz[a]
        if (!kept[a]) next
        mij <- mi_i[j]
        for (b in (a + 1L):length(idx_nz)) {
          k <- idx_nz[b]
          if (!kept[b]) next
          mik <- mi_i[k]
          # fetch MI(j, k)
          # to avoid reading entire row j, read a single value by slicing [j,k]
          mjk <- h5read(h5_path, dset, index = list(j, k))
          # DPI condition
          # if MI(j,k) >= min(MI(i,j), MI(i,k)) - tol, remove weaker edge (i,j) or (i,k)
          threshold <- min(mij, mik) - dpi_tolerance
          if (mjk >= threshold) {
            if (mij <= mik) {
              kept[a] <- FALSE
              # break inner loop since edge (i,j) removed
              break
            } else {
              kept[b] <- FALSE
            }
          }
        }
      }
    }
    # store kept neighbors with j > i to avoid duplicates
    nbrs <- idx_nz[kept]
    keep_i[[i]] <- nbrs[nbrs > i]
    rm(mi_i); gc(FALSE)
  }
  
  # Build sparse symmetric adjacency
  ii <- integer(0); jj <- integer(0)
  for (i in seq_len(p)) {
    if (length(keep_i[[i]]) > 0L) {
      ii <- c(ii, rep.int(i, length(keep_i[[i]])))
      jj <- c(jj, keep_i[[i]])
    }
  }
  if (length(ii) == 0L) {
    return(Matrix::Matrix(0, nrow = p, ncol = p, sparse = TRUE))
  }
  adj <- Matrix::sparseMatrix(i = c(ii, jj), j = c(jj, ii), x = 1, dims = c(p, p))
  adj
}

## ============================================================
## Build per-cell-type cache: file-backed MI + ARACNE + NST
## ============================================================
.ct_build_type_cache_fb <- function(lst_scrna,
                                    ct,
                                    cutoff,
                                    pcg,
                                    lrp,
                                    method_mi = "mm",
                                    ncores = NULL,
                                    echo = TRUE,
                                    cache_dir = tempdir(),
                                    block_size = 512,
                                    aracne_mode = c("stream", "full"),
                                    topk_per_node = Inf,
                                    dpi_tolerance = 0) {
  aracne_mode <- match.arg(aracne_mode)
  if (echo) message(sprintf("  [+] caching cell type: %s", ct))
  
  # Extract and filter
  mat_ct <- CytoTalk::extract_group(ct, lst_scrna)  # genes x cells (sparse)
  mat_filt <- CytoTalk::subset_rownames(.expr_filter(mat_ct, cutoff), pcg)
  if (nrow(mat_filt) == 0L) stop(sprintf("No genes passed cutoff/pcg for cell type '%s'.", ct))
  
  # Discretize (dense)
  dense_disc <- .ct_discretize_dense(mat_filt)
  p <- ncol(dense_disc)
  coln <- colnames(dense_disc); if (is.null(coln)) coln <- paste0("g", seq_len(p))
  
  # HDF5 path for MI
  h5_path <- file.path(cache_dir, sprintf("MI_%s.h5", gsub("[^A-Za-z0-9_]+", "_", ct)))
  # Compute MI to HDF5
  .compute_mi_to_hdf5(dense_disc, h5_path, dset = "MI", block_size = block_size, method_mi = method_mi, ncores = ncores)
  rm(dense_disc); gc(FALSE)
  
  # ARACNE from HDF5
  if (echo) message("     - ARACNE (", aracne_mode, ") ...")
  if (aracne_mode == "full") {
    # This will load the full MI (RAM spike). Use only if RAM permits.
    MI_full <- h5read(h5_path, "MI")
    MI_full <- .ct_zero_diag(MI_full)
    mat_intra <- Matrix::Matrix(parmigene::aracne.m(MI_full), sparse = TRUE)
    rm(MI_full); gc(FALSE)
  } else {
    # Streamed DPI pruning, no full MI in RAM
    mat_intra <- .aracne_from_hdf5(h5_path, dset = "MI", mode = "stream", topk_per_node = topk_per_node, dpi_tolerance = dpi_tolerance)
  }
  
  # Non-selftalk
  vec_nst <- CytoTalk::nonselftalk(mat_ct, lrp)
  
  list(
    mat_intra = mat_intra,
    vec_nst   = vec_nst,
    n_genes_kept = length(coln),
    mi_h5 = h5_path
  )
}

## ============================================================
## Main runner
## ============================================================
run_cytotalk_cached <- function(
    lst_scrna,
    cell_types,
    cutoff_a = 0.2,
    cutoff_b = 0.2,
    pcg = CytoTalk::pcg_human,
    lrp = CytoTalk::lrp_human,
    beta_max = 100,
    omega_min = 0.5,
    omega_max = 0.5,
    depth = 3,
    ntrial = 1000,
    cores = NULL,
    echo = TRUE,
    dir_out = NULL,
    ## NEW:
    cache_dir = "cytotalk_cache",     # directory for HDF5 MI files
    mi_block_size = 512,              # MI block size (tune for disk IO)
    aracne_mode = c("stream", "full"),
    topk_per_node = Inf,              # set finite (e.g. 1000) to speed DPI with large p
    dpi_tolerance = 0                 # DPI tolerance parameter
) {
  aracne_mode <- match.arg(aracne_mode)
  if (!is.null(dir_out) && !dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  
  # Step 1: PEM once
  if (echo) .ct_tick("Preprocessing (PEM & caches)...")
  mat_pem <- CytoTalk::pem(lst_scrna)
  
  # Cache env
  cache <- new.env(parent = emptyenv())
  
  # Decide cores globally
  if (is.null(cores) || cores < 1L) cores <- max(1L, parallel::detectCores() - 1L)
  
  # Build caches for each cell type to superset cutoff
  cutoff_superset <- max(cutoff_a, cutoff_b)
  for (ct in cell_types) {
    key <- paste0("ct:", ct)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      obj <- .ct_build_type_cache_fb(
        lst_scrna   = lst_scrna,
        ct          = ct,
        cutoff      = cutoff_superset,
        pcg         = pcg,
        lrp         = lrp,
        method_mi   = "mm",
        ncores      = cores,
        echo        = echo,
        cache_dir   = cache_dir,
        block_size  = mi_block_size,
        aracne_mode = aracne_mode,
        topk_per_node = topk_per_node,
        dpi_tolerance = dpi_tolerance
      )
      assign(key, obj, envir = cache)
      gc(FALSE)
    }
  }
  
  if (echo) .ct_tick("Ready (PEM + per-type MI/ARACNE cached).")
  
  # Pair runner
  run_pair <- function(type_a, type_b, beta_max_pair = beta_max) {
    if (echo) message(sprintf("Running CytoTalk (cached) for %s vs %s ...", type_a, type_b))
    obj_a <- get(paste0("ct:", type_a), envir = cache, inherits = FALSE)
    obj_b <- get(paste0("ct:", type_b), envir = cache, inherits = FALSE)
    mat_a <- CytoTalk::extract_group(type_a, lst_scrna)
    mat_b <- CytoTalk::extract_group(type_b, lst_scrna)
    
    if (echo) .ct_tick("Integrate network...")
    lst_net <- CytoTalk::integrate_network(
      obj_a$vec_nst, obj_b$vec_nst,
      obj_a$mat_intra, obj_b$mat_intra,
      type_a, type_b,
      mat_pem, mat_a, lrp
    )
    
    if (echo) .ct_tick("PCSF...")
    lst_pcst <- CytoTalk::run_pcst(lst_net, beta_max_pair, omega_min, omega_max)
    
    if (echo) .ct_tick("Determine best signaling network...")
    df_test <- CytoTalk::ks_test_pcst(lst_pcst)
    index <- order(as.numeric(df_test[, "pval"]))[1]
    beta <- df_test[index, "beta"]; omega <- df_test[index, "omega"]
    
    if (echo) .ct_tick("Generate network output...")
    df_net <- CytoTalk::extract_network(lst_net, lst_pcst, mat_pem, beta, omega)
    lst_path <- CytoTalk::extract_pathways(df_net, type_a, depth)
    lst_graph <- if (!is.null(lst_path)) lapply(lst_path, CytoTalk::graph_pathway) else NULL
    
    if (is.null(lst_path)) {
      if (echo) .ct_tick("NOTE: No pathways found, analysis skipped!")
      return(list(
        params = list(cell_type_a = type_a, cell_type_b = type_b,
                      beta_max = beta_max_pair, omega_min = omega_min,
                      omega_max = omega_max, depth = depth, ntrial = ntrial),
        pem = mat_pem,
        integrated_net = lst_net,
        pcst = list(occurances = lst_pcst, ks_test_pval = df_test, final_network = df_net),
        pathways = NULL
      ))
    }
    
    if (echo) .ct_tick("Analyze pathways...")
    lst_pval <- lapply(lst_path, CytoTalk::analyze_pathway, lst_net, type_a, type_b, beta, ntrial)
    
    nodes <- do.call(rbind, strsplit(names(lst_path), "--"))
    df_pval <- do.call(rbind, apply(nodes, 1, function(x) {
      t <- ifelse(endsWith(x, "_"), type_b, type_a)
      x <- gsub("_$", "", x)
      data.frame(
        ligand = x[1], receptor = x[2],
        ligand_type = t[1], receptor_type = t[2],
        stringsAsFactors = FALSE
      )
    }))
    df_pval <- cbind(df_pval, do.call(rbind, lst_pval))
    df_pval <- df_pval[order(as.numeric(df_pval$pval_potential)), , drop = FALSE]
    
    list(
      params = list(cell_type_a = type_a, cell_type_b = type_b,
                    beta_max = beta_max_pair, omega_min = omega_min,
                    omega_max = omega_max, depth = depth, ntrial = ntrial),
      pem = mat_pem,
      integrated_net = lst_net,
      pcst = list(occurances = lst_pcst, ks_test_pval = df_test, final_network = df_net),
      pathways = list(raw = lst_path, graphs = lst_graph, df_pval = df_pval)
    )
  }
  
  structure(
    list(
      run_pair = run_pair,
      cache    = cache,
      pem      = mat_pem,
      meta     = list(cache_dir = cache_dir, aracne_mode = aracne_mode, topk_per_node = topk_per_node)
    ),
    class = "cytotalk_cache_runner"
  )
}

runner <- run_cytotalk_cached(
  lst_scrna = lst_scrna,
  cell_types = levels(factor(lst_scrna$cell_types)),
  cutoff_a = 0.2,
  cutoff_b = 0.2,
  pcg = CytoTalk::pcg_human,
  lrp = CytoTalk::lrp_human,
  beta_max = 100,
  cores = NULL,
  echo = TRUE,
  cache_dir = "cytotalk_cache_h5",  # where MI files go
  mi_block_size = 512,              # tune if you have faster disks
  aracne_mode = "stream",           # << stays in low RAM
  topk_per_node = Inf,              # keep Inf to consider all neighbors for DPI
  dpi_tolerance = 0
)

# Example pair
res <- runner$run_pair("Astro", "EN_L2_3_IT", beta_max_pair = 100)
