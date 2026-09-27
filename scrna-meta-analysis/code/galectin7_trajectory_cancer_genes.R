# Galectin-7 versus cancer-cell programs and differentiation pseudotime.
#
# Uses the Seurat objects written by galectin7_oscc_scrna_seurat.R.
# Trajectory is fit inside tumor epithelial or malignant cells, rooted at the
# most basal cluster and aimed at the most differentiated cluster (slingshot).
# This is an expression-state ordering, not a measured clock.

PROJECT_DIR <- "F:/oral_cancer_galectin-7"
OUT_DIR <- file.path(PROJECT_DIR, "results", "scrna_galectin7", "trajectory")
RDS <- list(
  GSE181919 = file.path(PROJECT_DIR, "results", "scrna_galectin7", "GSE181919", "GSE181919_galectin.rds"),
  GSE103322 = file.path(PROJECT_DIR, "results", "scrna_galectin7", "GSE103322", "GSE103322_galectin.rds"),
  GSE172577 = file.path(PROJECT_DIR, "results", "scrna_galectin7", "GSE172577", "GSE172577_galectin.rds")
)

FOCUS <- "LGALS7"
MAX_TRAJ_CELLS <- 6000L

PROGRAMS <- list(
  basal = c("KRT15", "KRT5", "KRT14", "TP63", "COL17A1", "ITGA6"),
  differentiated = c("KRT13", "KRT4", "IVL", "SPRR1B", "SPINK5", "KRT1"),
  epithelial_adhesion = c("CDH1", "EPCAM", "DSG3", "PKP1", "JUP"),
  proliferation = c("MKI67", "TOP2A", "PCNA", "TYMS", "BIRC5", "CDK1"),
  emt = c("VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "CDH2", "FN1", "TWIST1"),
  invasion = c("MMP1", "MMP9", "MMP10", "LAMC2", "PTHLH", "CXCL8"),
  oncogenic_signaling = c("EGFR", "MYC", "CCND1", "CDKN2A", "STAT3", "JUN", "FOS")
)

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(ggplot2)
  library(dplyr)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(7)
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

seurat_v5 <- function() utils::packageVersion("Seurat") >= "5.0.0"

get_layer <- function(obj, layer) {
  if (seurat_v5()) GetAssayData(obj, assay = "RNA", layer = layer) else GetAssayData(obj, assay = "RNA", slot = layer)
}

gene_vec <- function(obj, gene) {
  mat <- get_layer(obj, "data")
  if (!gene %in% rownames(mat)) return(rep(NA_real_, ncol(obj)))
  as.numeric(mat[gene, ])
}

present_genes <- function(obj, genes) intersect(genes, rownames(obj))

tumor_epithelial_cells <- function(obj) {
  state <- as.character(obj$epithelial_state)
  tissue <- as.character(obj$tissue)
  ct <- as.character(obj$cell_type)
  use <- state %in% c("OSCC cancer cell", "LN cancer cell", "Tumor epithelial (marker-defined)")
  if (!any(use)) {
    use <- ct %in% c("Malignant", "Epithelial") & tissue == "primary_tumor"
  }
  colnames(obj)[use & !is.na(use)]
}

mean_program <- function(obj, genes) {
  genes <- present_genes(obj, genes)
  if (length(genes) == 0) return(rep(NA_real_, ncol(obj)))
  mat <- get_layer(obj, "data")
  as.numeric(colMeans(as.matrix(mat[genes, , drop = FALSE])))
}

spearman_row <- function(x, y, dataset, subset_name, feature, n_feature_genes = NA_integer_) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 30 || length(unique(x)) < 3 || length(unique(y)) < 3) return(NULL)
  test <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    dataset = dataset,
    subset = subset_name,
    feature = feature,
    n_cells = length(x),
    n_genes_in_feature = n_feature_genes,
    rho = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

correlate_programs <- function(obj, dataset, subset_name) {
  lg <- gene_vec(obj, FOCUS)
  rows <- list()
  for (nm in names(PROGRAMS)) {
    genes <- present_genes(obj, PROGRAMS[[nm]])
    if (length(genes) < 2) next
    score <- mean_program(obj, genes)
    rows[[length(rows) + 1L]] <- spearman_row(lg, score, dataset, subset_name, nm, length(genes))
    for (gene in genes) {
      rows[[length(rows) + 1L]] <- spearman_row(
        lg, gene_vec(obj, gene), dataset, subset_name, gene, 1L
      )
    }
  }
  other <- setdiff(intersect(c(
    "LGALS1", "LGALS2", "LGALS3", "LGALS4", "LGALS7B", "LGALS8", "LGALS9"
  ), rownames(obj)), FOCUS)
  for (gene in other) {
    rows[[length(rows) + 1L]] <- spearman_row(
      lg, gene_vec(obj, gene), dataset, subset_name, gene, 1L
    )
  }
  out <- bind_rows(rows)
  if (nrow(out) == 0) return(out)
  out$q_value <- stats::p.adjust(out$p_value, method = "BH")
  out
}

prepare_subset <- function(obj) {
  mode <- obj@misc$expr_mode
  if (identical(mode, "log2_tpm10")) {
    obj <- tryCatch(
      FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = 2000, verbose = FALSE),
      error = function(e) FindVariableFeatures(obj, selection.method = "mean.var.plot", nfeatures = 2000, verbose = FALSE)
    )
  } else {
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 20, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
  obj
}

choose_terminals <- function(obj) {
  avg <- AverageExpression(
    obj,
    features = unique(c(PROGRAMS$basal, PROGRAMS$differentiated)),
    group.by = "seurat_clusters",
    assays = "RNA"
  )[[1]]
  avg <- as.matrix(avg)
  colnames(avg) <- sub("^g([0-9].*)$", "\\1", colnames(avg))
  sizes <- table(obj$seurat_clusters)
  basal <- present_genes(obj, PROGRAMS$basal)
  diff <- present_genes(obj, PROGRAMS$differentiated)
  basal <- intersect(basal, rownames(avg))
  diff <- intersect(diff, rownames(avg))
  if (length(basal) == 0 || length(diff) == 0) {
    return(list(start = NA_character_, end = NA_character_))
  }
  score <- colMeans(avg[basal, , drop = FALSE]) - colMeans(avg[diff, , drop = FALSE])
  keep <- intersect(names(score), names(sizes)[sizes >= 25])
  score <- score[keep]
  score <- score[is.finite(score)]
  if (length(score) < 2) {
    return(list(start = NA_character_, end = NA_character_))
  }
  list(
    start = names(which.max(score)),
    end = names(which.min(score))
  )
}

run_slingshot <- function(obj, dataset) {
  terminals <- choose_terminals(obj)
  if (length(terminals$start) != 1L || is.na(terminals$start)) {
    stop("Not enough clusters for a basal-to-differentiated path.", call. = FALSE)
  }
  sce <- as.SingleCellExperiment(obj)
  rd_names <- reducedDimNames(sce)
  rd <- rd_names[grepl("PCA", rd_names, ignore.case = TRUE)]
  if (length(rd) == 0) stop("PCA reduction was not found.", call. = FALSE)
  rd <- rd[[1]]
  cl <- as.character(obj$seurat_clusters)
  fit <- tryCatch(
    slingshot(sce, clusterLabels = cl, reducedDim = rd, start.clus = terminals$start, end.clus = terminals$end),
    error = function(e) {
      message("Slingshot with an end cluster failed (", conditionMessage(e), "); retrying from the basal cluster only.")
      slingshot(sce, clusterLabels = cl, reducedDim = rd, start.clus = terminals$start)
    }
  )
  pt <- slingPseudotime(fit, na = TRUE)
  if (is.null(dim(pt))) pt <- matrix(pt, ncol = 1, dimnames = list(colnames(obj), "curve1"))
  # Prefer the curve that orders the most cells and spans basal to differentiated.
  n_ok <- apply(pt, 2, function(z) sum(is.finite(z)))
  curve <- names(which.max(n_ok))
  data.frame(
    cell = colnames(obj),
    dataset = dataset,
    sample_id = as.character(obj$sample_id),
    tissue = as.character(obj$tissue),
    cell_type = as.character(obj$cell_type),
    cluster = as.character(obj$seurat_clusters),
    pseudotime = as.numeric(pt[, curve]),
    curve = curve,
    n_curves = ncol(pt),
    start_cluster = terminals$start,
    end_cluster = terminals$end,
    LGALS7 = gene_vec(obj, FOCUS),
    basal = mean_program(obj, PROGRAMS$basal),
    differentiated = mean_program(obj, PROGRAMS$differentiated),
    stringsAsFactors = FALSE
  )
}

plot_pseudotime <- function(df, dataset) {
  d <- df[is.finite(df$pseudotime), , drop = FALSE]
  if (nrow(d) < 50) return(invisible(NULL))
  p <- ggplot(d, aes(pseudotime, LGALS7)) +
    geom_point(aes(color = differentiated - basal), alpha = 0.35, size = 0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "#B2182B", linewidth = 0.8) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "Diff - basal") +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(dataset, ": LGALS7 along differentiation pseudotime"),
      x = "Slingshot pseudotime (basal cluster as root)",
      y = "LGALS7 expression"
    )
  ggplot2::ggsave(file.path(OUT_DIR, paste0(dataset, "_LGALS7_pseudotime.pdf")), p, width = 7.2, height = 4.8)
  ggplot2::ggsave(file.path(OUT_DIR, paste0(dataset, "_LGALS7_pseudotime.png")), p, width = 7.2, height = 4.8, dpi = 160)
}

analyze_one <- function(dataset, rds) {
  message("\n======== ", dataset, " ========")
  obj <- readRDS(rds)
  cells <- tumor_epithelial_cells(obj)
  message("Tumor epithelial/malignant cells: ", length(cells))
  if (length(cells) < 150) stop("Too few tumor epithelial cells.", call. = FALSE)
  sub <- subset(obj, cells = cells)
  rm(obj)
  gc(verbose = FALSE)

  cor_all <- correlate_programs(sub, dataset, "all_tumor_epithelial")
  pos <- colnames(sub)[gene_vec(sub, FOCUS) > 0]
  cor_pos <- if (length(pos) >= 100) {
    correlate_programs(subset(sub, cells = pos), dataset, "LGALS7_detected")
  } else {
    cor_all[0, ]
  }

  traj_obj <- sub
  if (ncol(traj_obj) > MAX_TRAJ_CELLS) {
    set.seed(7)
    keep <- sample(colnames(traj_obj), MAX_TRAJ_CELLS)
    traj_obj <- subset(traj_obj, cells = keep)
    message("Downsampled trajectory to ", ncol(traj_obj), " cells")
  }
  traj_obj <- prepare_subset(traj_obj)
  pt <- run_slingshot(traj_obj, dataset)
  utils::write.csv(pt, file.path(OUT_DIR, paste0(dataset, "_pseudotime_cells.csv")), row.names = FALSE)
  plot_pseudotime(pt, dataset)

  pt_ok <- pt[is.finite(pt$pseudotime), , drop = FALSE]
  tests <- bind_rows(
    spearman_row(pt_ok$LGALS7, pt_ok$pseudotime, dataset, "trajectory_cells", "LGALS7_vs_pseudotime"),
    spearman_row(pt_ok$differentiated, pt_ok$pseudotime, dataset, "trajectory_cells", "differentiated_vs_pseudotime"),
    spearman_row(pt_ok$basal, pt_ok$pseudotime, dataset, "trajectory_cells", "basal_vs_pseudotime"),
    spearman_row(pt_ok$LGALS7, pt_ok$differentiated, dataset, "trajectory_cells", "LGALS7_vs_differentiated"),
    spearman_row(pt_ok$LGALS7, pt_ok$basal, dataset, "trajectory_cells", "LGALS7_vs_basal")
  )
  if (nrow(tests) > 0) tests$q_value <- stats::p.adjust(tests$p_value, method = "BH")

  rm(sub, traj_obj)
  gc(verbose = FALSE)
  list(cor = bind_rows(cor_all, cor_pos), tests = tests, n = nrow(pt_ok), start = pt$start_cluster[1], end = pt$end_cluster[1])
}

main_traj <- function() {
  cor_all <- list()
  tests_all <- list()
  notes <- character()
  for (dataset in names(RDS)) {
    result <- tryCatch(
      analyze_one(dataset, RDS[[dataset]]),
      error = function(e) {
        message("FAILED ", dataset, ": ", conditionMessage(e))
        structure(conditionMessage(e), class = "traj_error")
      }
    )
    if (inherits(result, "traj_error")) {
      notes <- c(notes, paste0(dataset, " FAILED: ", result))
      next
    }
    cor_all[[dataset]] <- result$cor
    tests_all[[dataset]] <- result$tests
    notes <- c(
      notes,
      paste0(
        dataset, ": trajectory n=", result$n,
        " start_cluster=", result$start, " end_cluster=", result$end
      )
    )
  }
  cor_df <- bind_rows(cor_all)
  test_df <- bind_rows(tests_all)
  if (is.null(cor_df)) cor_df <- data.frame()
  if (is.null(test_df)) test_df <- data.frame()
  if (nrow(cor_df) > 0) {
    utils::write.csv(cor_df, file.path(OUT_DIR, "LGALS7_cancer_program_correlations.csv"), row.names = FALSE)
  }
  if (nrow(test_df) > 0) {
    utils::write.csv(test_df, file.path(OUT_DIR, "LGALS7_pseudotime_tests.csv"), row.names = FALSE)
  }
  writeLines(notes, file.path(OUT_DIR, "RUN_NOTES.txt"))
  message("Done: ", OUT_DIR)
}

if (sys.nframe() == 0L) main_traj()
