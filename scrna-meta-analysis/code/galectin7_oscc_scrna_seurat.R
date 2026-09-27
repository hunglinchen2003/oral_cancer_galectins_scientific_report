# Galectin expression in oral squamous cell carcinoma scRNA-seq
#
# The files in oscc_raw_cpout/ and oscc_processed_data/ are imaging mass
# cytometry (IMC) from Nasiaho et al., Cell Reports Medicine 2026
# (doi:10.1016/j.xcrm.2026.102615). That assay measures 25 proteins
# (E-cadherin, KI67, PROX1, CD45, and so on). Galectin-7 is not in the panel,
# so those files cannot answer a galectin question.
#
# This script reanalyzes the three public scRNA-seq cohorts used in that
# paper, with the question shifted from proliferating lymphatic endothelial
# cells to the galectin family, especially galectin-7 (LGALS7):
#
#   GSE181919  Choi et al., Nat Commun 2023. 10x. Default subset is oral
#              cavity (subsite OC): 10 primary OSCC, 3 normal mucosa,
#              4 leukoplakia. Author cell labels include Malignant.cells
#              and Epithelial.cells. This is the cohort with a normal
#              oral epithelial comparator.
#   GSE103322  Puram et al., Cell 2017. Smart-seq2, 5,902 QC-filtered cells
#              from 18 oral-cavity HNSCC patients (primary tumor and lymph
#              node). Released values are already log2(TPM/10 + 1). Cancer
#              cells are the authors' CNV-based call.
#   GSE172577  Peng et al., Oral Oncol 2021. 10x, 6 OSCC tumors. No author
#              cell-type table; epithelial cells are marker-defined and are
#              not a CNV-confirmed malignant call.
#
# Expression values are not on the same scale across technologies. Compare
# LGALS7 within a dataset. Across datasets, compare percent of cells with
# a detected transcript, not the raw mean.
#
# Run from PowerShell:
#   & "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" analysis\galectin7_oscc_scrna_seurat.R
#
# Optional: set GALECTIN_DATASETS=GSE181919 to run one cohort.
# First run downloads GEO supplements (about 0.8 GB for all three cohorts).
# Oral-cavity GSE181919 is the preferred cohort for cancer vs normal.

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) {
    return(NA_character_)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

this_file <- script_path()
PROJECT_DIR <- if (!is.na(this_file)) {
  dirname(dirname(this_file))
} else {
  "F:/oral_cancer_galectin-7"
}

DATA_DIR <- file.path(PROJECT_DIR, "data", "geo_scrna")
OUT_DIR <- file.path(PROJECT_DIR, "results", "scrna_galectin7")

DEFAULT_DATASETS <- c("GSE181919", "GSE103322", "GSE172577")
env_datasets <- Sys.getenv("GALECTIN_DATASETS", unset = "")
DATASETS <- if (nzchar(env_datasets)) {
  trimws(strsplit(env_datasets, ",", fixed = TRUE)[[1]])
} else {
  DEFAULT_DATASETS
}

# GSE181919 also contains oropharynx and hypopharynx. Keep oral cavity.
CHOI_SUBSITE <- "OC"
INCLUDE_LEUKOPLAKIA <- TRUE
MIN_CELLS_PER_SAMPLE_TYPE <- 20
RUN_DECONTX <- FALSE
SAVE_RDS <- TRUE
N_PCS <- 30L
CLUSTER_RESOLUTION <- 0.5

GALECTIN_GENES <- c(
  "LGALS1", "LGALS2", "LGALS3", "LGALS4",
  "LGALS7", "LGALS7B",
  "LGALS8", "LGALS9", "LGALS9B", "LGALS9C",
  "LGALS12", "LGALS13", "LGALS14", "LGALS16"
)
FOCUS_GENE <- "LGALS7"

# LGALS7 is intentionally absent from the epithelial signature.
LINEAGE_MARKERS <- list(
  Epithelial = c("KRT5", "KRT14", "KRT6A", "KRT17", "EPCAM"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "FAP"),
  Endothelial = c("PECAM1", "CLDN5", "CDH5", "VWF"),
  Pericyte = c("RGS5", "ACTA2", "MCAM"),
  "T cell" = c("CD3D", "CD3E", "TRAC", "IL7R"),
  "B/Plasma" = c("MS4A1", "CD79A", "MZB1", "JCHAIN"),
  Myeloid = c("LYZ", "CD68", "LST1", "AIF1"),
  Mast = c("TPSAB1", "TPSB2", "KIT")
)

BASAL_GENES <- c("KRT15", "KRT5", "KRT14", "COL17A1")
DIFF_GENES <- c("KRT4", "KRT13", "IVL", "SPRR1B", "SPINK5")

CELLTYPE_ORDER <- c(
  "Malignant", "Epithelial", "Fibroblast", "Endothelial", "Pericyte",
  "T cell", "B cell", "B/Plasma", "Macrophage", "Dendritic", "Myeloid",
  "Mast", "Myocyte", "Unassigned"
)

CELLTYPE_COLORS <- c(
  Malignant = "#B2182B",
  Epithelial = "#EF8A62",
  Fibroblast = "#1B9E77",
  Endothelial = "#377EB8",
  Pericyte = "#984EA3",
  "T cell" = "#4DAF4A",
  "B cell" = "#FF7F00",
  "B/Plasma" = "#A6761D",
  Macrophage = "#E6AB02",
  Dendritic = "#F781BF",
  Myeloid = "#66C2A5",
  Mast = "#666666",
  Myocyte = "#A6CEE3",
  Unassigned = "#BDBDBD"
)

STATE_ORDER <- c(
  "OSCC cancer cell",
  "LN cancer cell",
  "Leukoplakia epithelial",
  "Normal epithelial",
  "Tumor epithelial (marker-defined)",
  "Other"
)

# 10x viability filter used for Peng, whose released matrices are Cell Ranger
# outputs rather than a finished cell-type table. Choi and Puram are already
# the author-filtered cell sets, so they are not filtered again.
PENG_MIN_FEATURES <- 200
PENG_MAX_FEATURES <- 7500
PENG_MAX_MT <- 20

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

required_packages <- c("Seurat", "Matrix", "data.table", "dplyr", "tidyr", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install these R packages before running:\n  install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(7)
data.table::setDTthreads(2)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}
options(future.globals.maxSize = 8 * 1024^3)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

seurat_v5 <- function() {
  utils::packageVersion("Seurat") >= "5.0.0"
}

get_layer <- function(obj, layer) {
  if (seurat_v5()) {
    GetAssayData(obj, assay = "RNA", layer = layer)
  } else {
    GetAssayData(obj, assay = "RNA", slot = layer)
  }
}

set_layer <- function(obj, layer, value) {
  if (seurat_v5()) {
    SetAssayData(obj, assay = "RNA", layer = layer, new.data = value)
  } else {
    SetAssayData(obj, assay = "RNA", slot = layer, new.data = value)
  }
}

join_if_needed <- function(obj) {
  tryCatch(JoinLayers(obj), error = function(e) obj)
}

gzip_complete <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(TRUE)
  }
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  ok <- TRUE
  repeat {
    chunk <- tryCatch(readBin(con, "raw", 8e6), error = function(e) {
      ok <<- FALSE
      raw(0)
    })
    if (!ok || length(chunk) == 0) {
      break
    }
  }
  ok
}

download_if_missing <- function(url, dest, min_bytes = 10000) {
  dest_ok <- function() {
    file.exists(dest) &&
      file.info(dest)$size >= min_bytes &&
      gzip_complete(dest)
  }
  if (dest_ok()) {
    message("Using existing file: ", dest)
    return(invisible(dest))
  }
  if (file.exists(dest)) {
    unlink(dest)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 3600)
  message("Downloading ", url)
  ok <- FALSE
  curl_bin <- Sys.which("curl")
  if (!nzchar(curl_bin)) {
    curl_bin <- Sys.which("curl.exe")
  }
  for (attempt in seq_len(3)) {
    status <- if (nzchar(curl_bin)) {
      system2(
        curl_bin,
        c(
          "-L", "--fail", "--retry", "2", "-C", "-",
          "--max-time", "3600", "-o", dest, url
        )
      )
    } else {
      tryCatch(
        utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
        error = function(e) {
          message("Download attempt ", attempt, " failed: ", conditionMessage(e))
          1
        }
      )
    }
    if (identical(as.integer(status), 0L) && dest_ok()) {
      ok <- TRUE
      break
    }
    message("Download attempt ", attempt, " did not produce a complete file.")
    if (file.exists(dest)) {
      unlink(dest)
    }
  }
  if (!ok) {
    stop("Download failed: ", url, call. = FALSE)
  }
  invisible(dest)
}

geo_suppl_url <- function(gse, filename) {
  stub <- paste0(substr(gse, 1, nchar(gse) - 3), "nnn")
  sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s",
    stub, gse, filename
  )
}

clean_name <- function(x) {
  x <- trimws(x)
  gsub("^'+|'+$", "", x)
}

to_sparse <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  Matrix::Matrix(x, sparse = TRUE)
}

sparse_from_numeric_dt <- function(dt, genes, chunk = 2500) {
  pieces <- list()
  n <- nrow(dt)
  for (start in seq(1, n, by = chunk)) {
    end <- min(start + chunk - 1L, n)
    block <- as.matrix(dt[start:end])
    storage.mode(block) <- "double"
    pieces[[length(pieces) + 1L]] <- Matrix::Matrix(block, sparse = TRUE)
    rm(block)
  }
  mat <- do.call(rbind, pieces)
  rownames(mat) <- genes
  mat
}

detect_expr_mode <- function(mat) {
  x <- if (inherits(mat, "sparseMatrix")) mat@x else as.numeric(mat)
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return("counts")
  }
  frac_int <- mean(abs(x - round(x)) < 1e-8)
  # Puram's released matrix is log2(TPM/10 + 1): small, non-integer values.
  if (max(x) < 50 && frac_int < 0.9) "log2_tpm10" else "counts"
}

log2tpm_to_tpm <- function(logmat) {
  tpm <- logmat
  if (inherits(tpm, "sparseMatrix")) {
    tpm@x <- (2 ^ tpm@x - 1) * 10
    tpm@x[tpm@x < 0] <- 0
    return(tpm)
  }
  tpm <- (2 ^ tpm - 1) * 10
  tpm[tpm < 0] <- 0
  tpm
}

colors_for <- function(levels) {
  cols <- CELLTYPE_COLORS[levels]
  cols[is.na(cols)] <- "#999999"
  stats::setNames(unname(cols), levels)
}

factor_celltype <- function(x) {
  present <- unique(as.character(x))
  ordered <- c(intersect(CELLTYPE_ORDER, present), setdiff(present, CELLTYPE_ORDER))
  factor(x, levels = ordered)
}

save_plot <- function(plot, path, width = 8, height = 6) {
  tryCatch(
    {
      ggplot2::ggsave(
        paste0(path, ".pdf"), plot,
        width = width, height = height, limitsize = FALSE
      )
      ggplot2::ggsave(
        paste0(path, ".png"), plot,
        width = width, height = height, dpi = 160, limitsize = FALSE
      )
      invisible(TRUE)
    },
    error = function(e) {
      message("Plot failed (", basename(path), "): ", conditionMessage(e))
      invisible(FALSE)
    }
  )
}

empty_meta <- function(cells, dataset) {
  data.frame(
    cell = cells,
    dataset = dataset,
    sample_id = NA_character_,
    patient_id = NA_character_,
    tissue = NA_character_,
    subsite = "oral_cavity",
    author_cell_type = NA_character_,
    cell_type = NA_character_,
    stringsAsFactors = FALSE,
    row.names = cells
  )
}

assign_epithelial_state <- function(meta) {
  ct <- as.character(meta$cell_type)
  tissue <- as.character(meta$tissue)
  state <- rep("Other", nrow(meta))
  state[ct == "Malignant" & tissue == "primary_tumor"] <- "OSCC cancer cell"
  state[ct == "Malignant" & tissue == "lymph_node"] <- "LN cancer cell"
  state[ct == "Epithelial" & tissue == "normal"] <- "Normal epithelial"
  state[ct == "Epithelial" & tissue == "leukoplakia"] <- "Leukoplakia epithelial"
  state[ct == "Epithelial" & tissue == "primary_tumor"] <- "Tumor epithelial (marker-defined)"
  meta$epithelial_state <- factor(state, levels = intersect(STATE_ORDER, unique(state)))
  meta
}

# ---------------------------------------------------------------------------
# GSE181919 (Choi) -- author labels, oral cavity by default
# ---------------------------------------------------------------------------

read_choi_metadata <- function(path) {
  meta <- data.table::fread(path, header = FALSE, skip = 1, quote = "")
  if (ncol(meta) < 9) {
    stop("Unexpected GSE181919 metadata layout.", call. = FALSE)
  }
  meta <- meta[, 1:9]
  data.table::setnames(
    meta,
    c(
      "barcode", "patient_id", "sample_id", "gender", "age",
      "tissue_code", "subsite_code", "hpv", "author_cell_type"
    )
  )
  meta[, barcode_dot := chartr("-", ".", barcode)]
  meta[, age := as.numeric(age)]
  meta[]
}

choi_celltype_label <- function(x) {
  map <- c(
    "Malignant.cells" = "Malignant",
    "Epithelial.cells" = "Epithelial",
    "T.cells" = "T cell",
    "B_Plasma.cells" = "B/Plasma",
    "Macrophages" = "Macrophage",
    "Endothelial.cells" = "Endothelial",
    "Fibroblasts" = "Fibroblast",
    "Dendritic.cells" = "Dendritic",
    "Mast.cells" = "Mast",
    "Myocytes" = "Myocyte"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

choi_tissue_label <- function(x) {
  map <- c(
    CA = "primary_tumor",
    NL = "normal",
    LP = "leukoplakia",
    LN = "lymph_node"
  )
  out <- unname(map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

load_gse181919 <- function() {
  meta_path <- download_if_missing(
    geo_suppl_url("GSE181919", "GSE181919_Barcode_metadata.txt.gz"),
    file.path(DATA_DIR, "GSE181919_Barcode_metadata.txt.gz"),
    min_bytes = 10000
  )
  count_path <- download_if_missing(
    geo_suppl_url("GSE181919", "GSE181919_UMI_counts.txt.gz"),
    file.path(DATA_DIR, "GSE181919_UMI_counts.txt.gz"),
    min_bytes = 1e6
  )

  meta <- read_choi_metadata(meta_path)
  if (!identical(CHOI_SUBSITE, "all")) {
    meta <- meta[subsite_code == CHOI_SUBSITE]
  }
  if (!isTRUE(INCLUDE_LEUKOPLAKIA)) {
    meta <- meta[tissue_code != "LP"]
  }
  if (nrow(meta) == 0) {
    stop("No GSE181919 cells left after the subsite filter.", call. = FALSE)
  }

  message("Reading GSE181919 count header")
  hdr <- clean_name(read_first_line(count_path))
  idx <- match(meta$barcode_dot, hdr)
  matched <- !is.na(idx)
  message(
    "Matched ", sum(matched), " / ", nrow(meta),
    " metadata barcodes to the count matrix"
  )
  if (mean(matched) < 0.9) {
    stop("Most GSE181919 barcodes did not match the count matrix.", call. = FALSE)
  }
  meta <- meta[matched]
  idx <- idx[matched]

  message("Reading GSE181919 UMI counts for ", nrow(meta), " cells")
  dt <- data.table::fread(
    count_path,
    header = FALSE,
    skip = 1,
    quote = "",
    select = c(1L, idx + 1L),
    nThread = 2,
    showProgress = TRUE
  )
  genes <- clean_name(dt[[1]])
  dt[, 1 := NULL]
  keep <- !duplicated(genes) & nzchar(genes)
  genes <- genes[keep]
  dt <- dt[keep]
  data.table::setnames(dt, meta$barcode_dot)
  mat <- sparse_from_numeric_dt(dt, genes)
  colnames(mat) <- meta$barcode_dot
  rm(dt)
  gc(verbose = FALSE)

  present <- intersect(GALECTIN_GENES, rownames(mat))
  message("Galectin genes in GSE181919: ", paste(present, collapse = ", "))
  if (!FOCUS_GENE %in% rownames(mat)) {
    stop("LGALS7 is absent from GSE181919.", call. = FALSE)
  }
  keep_genes <- Matrix::rowSums(mat > 0) >= 3
  keep_genes[intersect(GALECTIN_GENES, rownames(mat))] <- TRUE
  mat <- mat[keep_genes, , drop = FALSE]

  cell_meta <- empty_meta(colnames(mat), "GSE181919")
  cell_meta$sample_id <- meta$sample_id
  cell_meta$patient_id <- meta$patient_id
  cell_meta$tissue <- choi_tissue_label(meta$tissue_code)
  cell_meta$subsite <- ifelse(meta$subsite_code == "OC", "oral_cavity", meta$subsite_code)
  cell_meta$author_cell_type <- meta$author_cell_type
  cell_meta$cell_type <- choi_celltype_label(meta$author_cell_type)
  cell_meta$gender <- meta$gender
  cell_meta$age <- meta$age
  cell_meta$hpv <- meta$hpv

  list(counts = mat, meta = cell_meta, expr_mode = "counts")
}

read_first_line <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]
}

# ---------------------------------------------------------------------------
# GSE103322 (Puram) -- metadata rows sit above the genes
# ---------------------------------------------------------------------------

load_gse103322 <- function() {
  path <- download_if_missing(
    geo_suppl_url("GSE103322", "GSE103322_HNSCC_all_data.txt.gz"),
    file.path(DATA_DIR, "GSE103322_HNSCC_all_data.txt.gz"),
    min_bytes = 1e6
  )
  hdr <- clean_name(read_first_line(path))
  if (!nzchar(hdr[[1]])) {
    hdr <- hdr[-1]
  }
  message("Reading GSE103322 metadata rows")
  meta_dt <- data.table::fread(
    path,
    header = FALSE,
    skip = 1,
    nrows = 5,
    quote = "",
    colClasses = "character",
    nThread = 2,
    showProgress = FALSE
  )
  meta_names <- gsub("\\s+", " ", clean_name(meta_dt[[1]]))
  meta_dt[, 1 := NULL]
  if (ncol(meta_dt) != length(hdr)) {
    stop(
      "GSE103322 metadata columns (", ncol(meta_dt),
      ") do not match cell names (", length(hdr), ").",
      call. = FALSE
    )
  }
  data.table::setnames(meta_dt, as.character(hdr))

  message("Reading GSE103322 expression matrix")
  expr_dt <- data.table::fread(
    path,
    header = FALSE,
    skip = 6,
    quote = "",
    nThread = 2,
    showProgress = TRUE
  )
  genes <- clean_name(expr_dt[[1]])
  expr_dt[, 1 := NULL]
  if (ncol(expr_dt) != length(hdr)) {
    stop("GSE103322 expression columns do not match the header.", call. = FALSE)
  }
  keep <- !duplicated(genes) & nzchar(genes)
  genes <- genes[keep]
  expr_dt <- expr_dt[keep]
  if ("Lymph node" %in% genes) {
    stop("GSE103322 parser left a metadata row in the gene matrix.", call. = FALSE)
  }
  if (!FOCUS_GENE %in% genes) {
    stop("LGALS7 is absent from GSE103322.", call. = FALSE)
  }
  data.table::setnames(expr_dt, hdr)
  logmat <- sparse_from_numeric_dt(expr_dt, genes)
  colnames(logmat) <- hdr
  rm(expr_dt)
  gc(verbose = FALSE)

  mode <- detect_expr_mode(logmat)
  message("GSE103322 expression mode: ", mode)
  if (mode != "log2_tpm10") {
    warning(
      "GSE103322 values were expected to be log2(TPM/10+1). ",
      "The detector called them ", mode, ".",
      call. = FALSE
    )
  }

  meta_row <- function(label) {
    hit <- which(meta_names == label)
    if (length(hit) != 1) {
      stop(
        "Missing Puram metadata row: ", label,
        "\nRows found: ", paste(meta_names, collapse = " | "),
        call. = FALSE
      )
    }
    trimws(as.character(unlist(meta_dt[hit], use.names = FALSE)))
  }
  # The released row name contains a double space.
  cancer <- as.integer(meta_row("classified as cancer cell")) == 1L
  lymph_node <- as.integer(meta_row("Lymph node")) == 1L
  raw_type <- meta_row("non-cancer cell type")
  raw_type <- sub("^-+", "", trimws(raw_type))
  raw_type[raw_type %in% c("0", "0.0", "", "NA")] <- NA_character_
  raw_type[tolower(raw_type) == "myocyte"] <- "Myocyte"

  cell_meta <- empty_meta(hdr, "GSE103322")
  cell_meta$author_cell_type <- raw_type
  cell_meta$author_cell_type[cancer] <- "Malignant"
  cell_meta$cell_type <- cell_meta$author_cell_type
  cell_meta$cell_type[is.na(cell_meta$cell_type)] <- "Unassigned"
  cell_meta$tissue <- ifelse(lymph_node, "lymph_node", "primary_tumor")
  cell_meta$patient_id <- sub("^(HNSCC[0-9]+|HN[0-9]+).*", "\\1", hdr)
  cell_meta$sample_id <- paste(cell_meta$patient_id, cell_meta$tissue, sep = "_")
  cell_meta$maxima_rt <- as.integer(meta_row("processed by Maxima enzyme"))

  list(counts = logmat, meta = cell_meta, expr_mode = mode)
}

# ---------------------------------------------------------------------------
# GSE172577 (Peng) -- 10x matrices inside RAW.tar
# ---------------------------------------------------------------------------

load_gse172577 <- function() {
  tar_path <- download_if_missing(
    geo_suppl_url("GSE172577", "GSE172577_RAW.tar"),
    file.path(DATA_DIR, "GSE172577_RAW.tar"),
    min_bytes = 1e8
  )
  exdir <- file.path(DATA_DIR, "GSE172577_RAW")
  if (!dir.exists(exdir) || length(list.files(exdir, recursive = TRUE)) < 6) {
    message("Extracting ", tar_path)
    untar(tar_path, exdir = exdir)
  }
  mtx_files <- list.files(exdir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE, recursive = TRUE)
  if (length(mtx_files) == 0) {
    stop("No 10x matrix files found in GSE172577_RAW.tar.", call. = FALSE)
  }

  mats <- list()
  metas <- list()
  for (mtx in mtx_files) {
    bn <- basename(mtx)
    parsed <- regexec("^(GSM[0-9]+)_(.+)_matrix\\.mtx\\.gz$", bn)
    hit <- regmatches(bn, parsed)[[1]]
    if (length(hit) != 3) {
      stop("Cannot parse sample name from ", bn, call. = FALSE)
    }
    sample_id <- hit[[3]]
    prefix <- paste0(hit[[2]], "_", sample_id)
    folder <- dirname(mtx)
    barcodes <- file.path(folder, paste0(prefix, "_barcodes.tsv.gz"))
    features <- file.path(folder, paste0(prefix, "_features.tsv.gz"))
    if (!file.exists(features)) {
      features <- file.path(folder, paste0(prefix, "_genes.tsv.gz"))
    }
    message("Reading ", sample_id)
    mat <- ReadMtx(
      mtx = mtx,
      cells = barcodes,
      features = features,
      feature.column = 2
    )
    rownames(mat) <- make.unique(rownames(mat))
    colnames(mat) <- paste(sample_id, colnames(mat), sep = "_")
    meta <- empty_meta(colnames(mat), "GSE172577")
    meta$sample_id <- sample_id
    meta$patient_id <- sample_id
    meta$tissue <- "primary_tumor"
    meta$subsite <- "oral_cavity"
    mats[[sample_id]] <- mat
    metas[[sample_id]] <- meta
  }

  genes <- Reduce(intersect, lapply(mats, rownames))
  mats <- lapply(mats, function(mat) mat[genes, , drop = FALSE])
  counts <- do.call(cbind, mats)
  meta <- do.call(rbind, metas)
  meta <- meta[colnames(counts), , drop = FALSE]
  if (!FOCUS_GENE %in% rownames(counts)) {
    stop("LGALS7 is absent from GSE172577.", call. = FALSE)
  }
  list(counts = counts, meta = meta, expr_mode = "counts", apply_qc = TRUE)
}

# ---------------------------------------------------------------------------
# Seurat workflow
# ---------------------------------------------------------------------------

repair_features <- function(mat) {
  # Seurat rejects underscores in feature names and rewrites them to dashes.
  rownames(mat) <- make.unique(gsub("_", "-", rownames(mat), fixed = TRUE))
  mat
}

align_meta <- function(meta, cells) {
  if (nrow(meta) != length(cells)) {
    stop("Metadata rows (", nrow(meta), ") do not match cells (", length(cells), ").", call. = FALSE)
  }
  if (identical(rownames(meta), cells)) {
    return(meta)
  }
  if (setequal(rownames(meta), cells)) {
    return(meta[cells, , drop = FALSE])
  }
  message("Cell metadata was not named by barcode; keeping the matrix column order.")
  rownames(meta) <- cells
  meta
}

build_seurat <- function(loaded, dataset) {
  counts <- repair_features(loaded$counts)
  meta <- align_meta(loaded$meta, colnames(counts))
  mode <- loaded$expr_mode
  if (mode == "log2_tpm10") {
    tpm <- log2tpm_to_tpm(counts)
    logmat <- counts
    obj <- CreateSeuratObject(
      counts = tpm,
      meta.data = meta,
      project = dataset,
      min.cells = 0,
      min.features = 0
    )
    obj <- join_if_needed(obj)
    if (nrow(logmat) != nrow(obj) || ncol(logmat) != ncol(obj)) {
      stop("The log expression matrix does not match the Seurat object.", call. = FALSE)
    }
    rownames(logmat) <- rownames(obj)
    colnames(logmat) <- colnames(obj)
    obj <- set_layer(obj, "data", logmat)
  } else {
    obj <- CreateSeuratObject(
      counts = counts,
      meta.data = meta,
      project = dataset,
      min.cells = 0,
      min.features = 0
    )
    obj <- join_if_needed(obj)
  }
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  if (all(is.na(obj$sample_id))) {
    parsed <- sub("_.*$", "", colnames(obj))
    if (dataset == "GSE172577" && any(grepl("^SYSMH", parsed))) {
      message("Restoring Peng sample labels from cell names.")
      obj$dataset <- dataset
      obj$sample_id <- parsed
      obj$patient_id <- parsed
      obj$tissue <- "primary_tumor"
      obj$subsite <- "oral_cavity"
    }
  }
  obj@misc$expression_scale <- if (mode == "log2_tpm10") {
    "log2(TPM/10+1)"
  } else {
    "log1p of counts scaled to 10,000"
  }
  obj@misc$expr_mode <- mode
  obj
}

apply_peng_qc <- function(obj) {
  before <- ncol(obj)
  obj <- subset(
    obj,
    subset = nFeature_RNA >= PENG_MIN_FEATURES &
      nFeature_RNA <= PENG_MAX_FEATURES &
      percent.mt < PENG_MAX_MT
  )
  message(
    "Peng QC kept ", ncol(obj), " / ", before,
    " cells (nFeature ", PENG_MIN_FEATURES, "-", PENG_MAX_FEATURES,
    ", percent.mt < ", PENG_MAX_MT, ")"
  )
  obj
}

reduce_seurat <- function(obj, batch_var = "sample_id") {
  mode <- obj@misc$expr_mode
  if (identical(mode, "counts")) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  } else {
    obj <- tryCatch(
      FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = 2000, verbose = FALSE),
      error = function(e) {
        message("dispersion variable-feature selection failed; using mean.var.plot")
        FindVariableFeatures(obj, selection.method = "mean.var.plot", nfeatures = 2000, verbose = FALSE)
      }
    )
  }
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = max(40L, N_PCS), verbose = FALSE)
  reduction <- "pca"
  batches <- unique(as.character(obj[[batch_var, drop = TRUE]]))
  batches <- batches[!is.na(batches)]
  if (length(batches) > 1 && requireNamespace("harmony", quietly = TRUE)) {
    obj <- tryCatch(
      {
        harmony::RunHarmony(
          obj,
          group.by.vars = batch_var,
          reduction.use = "pca",
          dims.use = 1:N_PCS
        )
      },
      error = function(e) {
        message("Harmony was not applied: ", conditionMessage(e))
        obj
      }
    )
    if ("harmony" %in% Reductions(obj)) {
      reduction <- "harmony"
    }
  } else if (length(batches) > 1) {
    message("harmony is not installed. UMAP uses PCA and can mix patients.")
  }
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:N_PCS, verbose = FALSE)
  obj <- FindClusters(obj, resolution = CLUSTER_RESOLUTION, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = reduction, dims = 1:N_PCS, verbose = FALSE)
  obj@misc$reduction <- reduction
  obj
}

add_marker_lineage <- function(obj) {
  available <- lapply(LINEAGE_MARKERS, function(genes) intersect(genes, rownames(obj)))
  available <- available[vapply(available, length, integer(1)) >= 2]
  if (length(available) == 0) {
    obj$marker_lineage <- NA_character_
    return(obj)
  }
  obj <- AddModuleScore(obj, features = available, name = "lineageScore", search = FALSE)
  score_cols <- paste0("lineageScore", seq_along(available))
  scores <- obj@meta.data[, score_cols, drop = FALSE]
  best <- apply(scores, 1, which.max)
  best_score <- apply(scores, 1, max)
  label <- names(available)[best]
  label[best_score < 0] <- "Unassigned"
  obj$marker_lineage <- label
  drop_cols <- grep("^lineageScore[0-9]+$", colnames(obj@meta.data), value = TRUE)
  obj@meta.data <- obj@meta.data[, setdiff(colnames(obj@meta.data), drop_cols), drop = FALSE]
  obj
}

add_keratinocyte_layer <- function(obj) {
  obj$keratinocyte_layer <- NA_character_
  use <- as.character(obj$epithelial_state) != "Other"
  if (sum(use) < 50) {
    return(obj)
  }
  if (sum(BASAL_GENES %in% rownames(obj)) < 2 || sum(DIFF_GENES %in% rownames(obj)) < 2) {
    return(obj)
  }
  sub <- subset(obj, cells = colnames(obj)[use])
  sub <- AddModuleScore(
    sub,
    features = list(basal = BASAL_GENES, differentiated = DIFF_GENES),
    name = "layerScore",
    search = TRUE
  )
  layer <- ifelse(sub$layerScore2 > sub$layerScore1, "Differentiated", "Basal")
  obj$keratinocyte_layer[colnames(sub)] <- layer
  obj
}

# ---------------------------------------------------------------------------
# Galectin summaries and plots
# ---------------------------------------------------------------------------

gene_vector <- function(obj, gene, layer) {
  values <- get_layer(obj, layer)
  if (!gene %in% rownames(values)) {
    return(rep(NA_real_, ncol(obj)))
  }
  as.numeric(values[gene, ])
}

cell_level_table <- function(obj) {
  meta <- obj@meta.data
  meta$cell <- rownames(meta)
  meta$LGALS7 <- gene_vector(obj, FOCUS_GENE, "data")
  meta$LGALS7_detected <- gene_vector(obj, FOCUS_GENE, "counts") > 0
  keep <- intersect(
    c(
      "cell", "dataset", "sample_id", "patient_id", "tissue", "subsite",
      "author_cell_type", "cell_type", "marker_lineage", "epithelial_state",
      "keratinocyte_layer", "nCount_RNA", "nFeature_RNA", "percent.mt",
      "LGALS7", "LGALS7_detected", "seurat_clusters"
    ),
    colnames(meta)
  )
  meta[, keep, drop = FALSE]
}

summarise_galectins <- function(obj) {
  genes <- intersect(GALECTIN_GENES, rownames(obj))
  meta <- obj@meta.data
  pieces <- lapply(genes, function(gene) {
    expr <- gene_vector(obj, gene, "data")
    detected <- gene_vector(obj, gene, "counts") > 0
    meta$expr <- expr
    meta$detected <- detected
    meta %>%
      group_by(dataset, sample_id, patient_id, tissue, subsite, cell_type, epithelial_state) %>%
      summarise(
        n_cells = dplyr::n(),
        n_detected = sum(detected, na.rm = TRUE),
        pct_detected = 100 * mean(detected, na.rm = TRUE),
        mean_expr = mean(expr, na.rm = TRUE),
        median_expr = median(expr, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(gene = gene, expression_scale = obj@misc$expression_scale)
  })
  bind_rows(pieces)
}

wilcox_row <- function(a, b, dataset, comparison, metric) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) {
    return(NULL)
  }
  test <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE))
  data.frame(
    dataset = dataset,
    comparison = comparison,
    metric = metric,
    n_a = length(a),
    n_b = length(b),
    median_a = median(a),
    median_b = median(b),
    mean_a = mean(a),
    mean_b = mean(b),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

sample_level_tests <- function(summary_tbl) {
  focus <- summary_tbl %>%
    filter(gene == FOCUS_GENE, n_cells >= MIN_CELLS_PER_SAMPLE_TYPE)
  if (nrow(focus) == 0) {
    return(focus[0, ])
  }
  dataset <- unique(as.character(focus$dataset))[[1]]
  tests <- list()

  state_metric <- function(state, metric) {
    focus %>%
      filter(as.character(epithelial_state) == state) %>%
      group_by(sample_id) %>%
      summarise(value = mean(.data[[metric]]), .groups = "drop") %>%
      pull(value)
  }
  comparisons <- list(
    c("OSCC cancer cell", "Normal epithelial"),
    c("OSCC cancer cell", "Leukoplakia epithelial"),
    c("OSCC cancer cell", "LN cancer cell")
  )
  for (pair in comparisons) {
    for (metric in c("mean_expr", "pct_detected")) {
      tests[[length(tests) + 1L]] <- wilcox_row(
        state_metric(pair[[1]], metric),
        state_metric(pair[[2]], metric),
        dataset,
        paste(pair[[1]], "vs", pair[[2]]),
        metric
      )
    }
  }

  tumor <- focus %>% filter(tissue == "primary_tumor")
  ref <- if ("Malignant" %in% tumor$cell_type) "Malignant" else if ("Epithelial" %in% tumor$cell_type) "Epithelial" else NA_character_
  others <- setdiff(unique(as.character(tumor$cell_type)), c(ref, "Unassigned"))
  if (!is.na(ref)) {
    for (other in others) {
      for (metric in c("mean_expr", "pct_detected")) {
        wide <- tumor %>%
          filter(cell_type %in% c(ref, other)) %>%
          select(sample_id, cell_type, all_of(metric)) %>%
          pivot_wider(names_from = cell_type, values_from = all_of(metric))
        if (!all(c(ref, other) %in% colnames(wide))) {
          next
        }
        paired <- wide %>% filter(is.finite(.data[[ref]]), is.finite(.data[[other]]))
        if (nrow(paired) < 3) {
          next
        }
        test <- suppressWarnings(stats::wilcox.test(paired[[ref]], paired[[other]], paired = TRUE, exact = FALSE))
        tests[[length(tests) + 1L]] <- data.frame(
          dataset = dataset,
          comparison = paste0("paired primary tumor: ", ref, " vs ", other),
          metric = metric,
          n_a = nrow(paired),
          n_b = nrow(paired),
          median_a = median(paired[[ref]]),
          median_b = median(paired[[other]]),
          mean_a = mean(paired[[ref]]),
          mean_b = mean(paired[[other]]),
          p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- bind_rows(tests)
  if (is.null(out) || !is.data.frame(out) || nrow(out) == 0) {
    return(data.frame())
  }
  out
}

plot_qc <- function(obj, fig_dir) {
  qc <- obj@meta.data
  qc$cell <- rownames(qc)
  long <- qc %>%
    select(sample_id, nFeature_RNA, nCount_RNA, percent.mt) %>%
    pivot_longer(-sample_id, names_to = "metric", values_to = "value")
  p <- ggplot(long, aes(sample_id, value)) +
    geom_violin(scale = "width", fill = "#d9d9d9") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = NULL, title = "QC metrics")
  save_plot(p, file.path(fig_dir, "qc_metrics"), width = 8, height = 9)
}

plot_umaps <- function(obj, fig_dir) {
  cols <- colors_for(levels(factor_celltype(obj$cell_type)))
  p_type <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE, cols = cols) +
    ggtitle(paste0(obj@project.name, " cell type"))
  save_plot(p_type, file.path(fig_dir, "umap_cell_type"), width = 9, height = 6)
  if (length(unique(obj$tissue)) > 1) {
    p_tissue <- DimPlot(obj, group.by = "tissue") +
      ggtitle(paste0(obj@project.name, " tissue"))
    save_plot(p_tissue, file.path(fig_dir, "umap_tissue"), width = 8, height = 6)
  }
  p_gene <- FeaturePlot(obj, features = FOCUS_GENE, order = TRUE) +
    ggtitle(paste0(FOCUS_GENE, "  [", obj@misc$expression_scale, "]"))
  save_plot(p_gene, file.path(fig_dir, "umap_LGALS7"), width = 7, height = 6)
}

plot_expression <- function(obj, summary_tbl, fig_dir) {
  obj$cell_type <- factor_celltype(obj$cell_type)
  cols <- colors_for(levels(obj$cell_type))
  p_vln <- VlnPlot(obj, features = FOCUS_GENE, group.by = "cell_type", pt.size = 0, cols = cols) +
    ggtitle(paste0(FOCUS_GENE, " by cell type")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_vln, file.path(fig_dir, "violin_LGALS7_by_cell_type"), width = 9, height = 5)

  genes <- intersect(GALECTIN_GENES, rownames(obj))
  p_dot <- DotPlot(obj, features = genes, group.by = "cell_type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Galectin family")
  save_plot(p_dot, file.path(fig_dir, "dotplot_galectins_by_cell_type"), width = 10, height = 6)

  focus_cells <- summary_tbl %>%
    filter(gene == FOCUS_GENE) %>%
    group_by(cell_type, tissue) %>%
    summarise(
      pct_detected = weighted.mean(pct_detected, n_cells),
      n_cells = sum(n_cells),
      .groups = "drop"
    ) %>%
    filter(n_cells >= 30)
  if (nrow(focus_cells) > 0) {
    focus_cells$cell_type <- factor_celltype(focus_cells$cell_type)
    p_pct <- ggplot(focus_cells, aes(cell_type, pct_detected, fill = cell_type)) +
      geom_col() +
      scale_fill_manual(values = colors_for(levels(focus_cells$cell_type))) +
      facet_wrap(~tissue, scales = "free_x") +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(
        x = NULL, y = paste0("% cells with ", FOCUS_GENE, " detected"),
        title = "Cell-weighted detection rate"
      )
    save_plot(p_pct, file.path(fig_dir, "LGALS7_detection_by_cell_type"), width = 10, height = 5)
  }

  per_sample <- summary_tbl %>%
    filter(gene == FOCUS_GENE, n_cells >= MIN_CELLS_PER_SAMPLE_TYPE)
  if (nrow(per_sample) > 0) {
    per_sample$cell_type <- factor_celltype(per_sample$cell_type)
    p_sample <- ggplot(per_sample, aes(cell_type, pct_detected, color = tissue)) +
      geom_boxplot(outlier.shape = NA, color = "grey40") +
      geom_jitter(width = 0.15, size = 1.8, alpha = 0.85) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = NULL, y = paste0("% ", FOCUS_GENE, "+ cells"),
        title = "One point per sample (groups with <20 cells omitted)"
      )
    save_plot(p_sample, file.path(fig_dir, "LGALS7_detection_per_sample"), width = 10, height = 5)
  }

  positive <- cell_level_table(obj) %>%
    filter(LGALS7_detected %in% TRUE, tissue == "primary_tumor")
  if (nrow(positive) > 0) {
    comp <- positive %>%
      count(sample_id, cell_type, name = "n_pos") %>%
      group_by(sample_id) %>%
      mutate(fraction = n_pos / sum(n_pos)) %>%
      ungroup()
    comp$cell_type <- factor_celltype(comp$cell_type)
    p_comp <- ggplot(comp, aes(sample_id, fraction, fill = cell_type)) +
      geom_col() +
      scale_fill_manual(values = colors_for(levels(comp$cell_type))) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        x = NULL, y = paste0("Fraction of ", FOCUS_GENE, "+ cells"),
        title = paste0("Which cells account for ", FOCUS_GENE, "+ transcripts in primary tumors")
      )
    save_plot(p_comp, file.path(fig_dir, "LGALS7_positive_composition"), width = 9, height = 5)
  }

  states <- summary_tbl %>%
    filter(gene == FOCUS_GENE, as.character(epithelial_state) != "Other", n_cells >= 10)
  if (nrow(states) > 0) {
    p_state <- ggplot(states, aes(epithelial_state, mean_expr, color = epithelial_state)) +
      geom_boxplot(outlier.shape = NA, color = "grey40") +
      geom_jitter(width = 0.12, size = 2) +
      theme_classic(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none") +
      labs(
        x = NULL,
        y = paste0("Mean ", FOCUS_GENE, " [", unique(states$expression_scale)[[1]], "]"),
        title = "One point per sample"
      )
    save_plot(p_state, file.path(fig_dir, "LGALS7_epithelial_states"), width = 8, height = 5)

    epi_cells <- colnames(obj)[as.character(obj$epithelial_state) != "Other"]
    if (length(epi_cells) >= 50) {
      sub <- subset(obj, cells = epi_cells)
      p_epi <- DotPlot(sub, features = genes, group.by = "epithelial_state") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Galectins in epithelial and cancer cells")
      save_plot(p_epi, file.path(fig_dir, "dotplot_galectins_epithelial_states"), width = 10, height = 5)
    }
  }

  if ("keratinocyte_layer" %in% colnames(obj@meta.data)) {
    layer_tbl <- cell_level_table(obj) %>%
      filter(!is.na(keratinocyte_layer), as.character(epithelial_state) != "Other") %>%
      group_by(epithelial_state, keratinocyte_layer, sample_id) %>%
      summarise(
        n_cells = dplyr::n(),
        mean_expr = mean(LGALS7),
        pct_detected = 100 * mean(LGALS7_detected),
        .groups = "drop"
      ) %>%
      filter(n_cells >= 10)
    if (nrow(layer_tbl) > 0) {
      p_layer <- ggplot(layer_tbl, aes(keratinocyte_layer, mean_expr, color = keratinocyte_layer)) +
        geom_boxplot(outlier.shape = NA, color = "grey40") +
        geom_jitter(width = 0.1, size = 1.6) +
        facet_wrap(~epithelial_state, scales = "free_y") +
        theme_classic(base_size = 11) +
        theme(legend.position = "none") +
        labs(x = NULL, y = paste0("Mean ", FOCUS_GENE), title = "Basal vs differentiated keratinocyte program")
      save_plot(p_layer, file.path(fig_dir, "LGALS7_keratinocyte_layer"), width = 9, height = 5)
    }
  }
}

write_dataset_summary <- function(obj, summary_tbl, tests, path) {
  focus <- summary_tbl %>% filter(gene == FOCUS_GENE)
  pooled <- focus %>%
    group_by(tissue, cell_type) %>%
    summarise(
      pct_detected = round(weighted.mean(pct_detected, n_cells), 1),
      mean_expr = round(weighted.mean(mean_expr, n_cells), 3),
      n_cells = sum(n_cells),
      .groups = "drop"
    ) %>%
    arrange(tissue, desc(pct_detected))
  lines <- c(
    paste0("Dataset: ", unique(obj$dataset)),
    paste0("Cells: ", ncol(obj)),
    paste0("Expression scale: ", obj@misc$expression_scale),
    paste0("Reduction: ", obj@misc$reduction),
    paste0("LGALS7 definition of detected: counts slot > 0"),
    "",
    "Cell-weighted LGALS7 by tissue and cell type:",
    capture.output(print(as.data.frame(pooled), row.names = FALSE)),
    "",
    "Sample-level tests use one value per sample, not one value per cell.",
    "Paired tests require at least 20 cells of each type in a sample and at least 3 samples."
  )
  if (nrow(tests) > 0) {
    lines <- c(lines, "", capture.output(print(tests, row.names = FALSE)))
  } else {
    lines <- c(lines, "", "No sample-level test met the sample-size rule.")
  }
  writeLines(lines, path)
}

analyze_dataset <- function(dataset) {
  message("\n======== ", dataset, " ========")
  loaded <- switch(
    dataset,
    GSE181919 = load_gse181919(),
    GSE103322 = load_gse103322(),
    GSE172577 = load_gse172577(),
    stop("Unknown dataset: ", dataset, call. = FALSE)
  )
  obj <- build_seurat(loaded, dataset)
  rm(loaded)
  gc(verbose = FALSE)

  out_dir <- file.path(OUT_DIR, dataset)
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  plot_qc(obj, fig_dir)

  if (identical(dataset, "GSE172577")) {
    obj <- apply_peng_qc(obj)
  }
  obj <- reduce_seurat(obj, batch_var = "sample_id")
  obj <- add_marker_lineage(obj)
  if (all(is.na(obj$cell_type))) {
    obj$cell_type <- obj$marker_lineage
    obj$author_cell_type <- "marker_module_score"
  }
  obj$cell_type <- factor_celltype(obj$cell_type)
  obj@meta.data <- assign_epithelial_state(obj@meta.data)
  obj <- add_keratinocyte_layer(obj)

  if (!is.na(obj$author_cell_type[[1]]) && !all(obj$author_cell_type == "marker_module_score")) {
    agree <- as.data.frame(table(author = obj$cell_type, marker = obj$marker_lineage))
    utils::write.csv(agree, file.path(out_dir, "author_vs_marker_counts.csv"), row.names = FALSE)
  }

  cells <- cell_level_table(obj)
  utils::write.csv(cells, file.path(out_dir, "cell_metadata_LGALS7.csv"), row.names = FALSE)
  summary_tbl <- summarise_galectins(obj)
  utils::write.csv(summary_tbl, file.path(out_dir, "galectin_by_sample_celltype.csv"), row.names = FALSE)
  tests <- sample_level_tests(summary_tbl)
  if (is.data.frame(tests) && nrow(tests) > 0) {
    utils::write.csv(tests, file.path(out_dir, "sample_level_tests.csv"), row.names = FALSE)
  }
  plot_umaps(obj, fig_dir)
  plot_expression(obj, summary_tbl, fig_dir)
  write_dataset_summary(obj, summary_tbl, tests, file.path(out_dir, "SUMMARY.txt"))
  if (isTRUE(SAVE_RDS)) {
    saveRDS(obj, file.path(out_dir, paste0(dataset, "_galectin.rds")))
  }
  message("Finished ", dataset, " -> ", out_dir)
  invisible(obj)
}

combine_datasets <- function() {
  files <- list.files(
    OUT_DIR,
    pattern = "^galectin_by_sample_celltype\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0) {
    return(invisible(NULL))
  }
  combined <- bind_rows(lapply(files, utils::read.csv, stringsAsFactors = FALSE))
  focus <- combined %>%
    filter(gene == FOCUS_GENE) %>%
    group_by(dataset, tissue, cell_type, expression_scale) %>%
    summarise(
      n_samples = dplyr::n_distinct(sample_id),
      pct_detected = weighted.mean(pct_detected, n_cells),
      n_cells = sum(n_cells),
      .groups = "drop"
    )
  out <- file.path(OUT_DIR, "combined")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(focus, file.path(out, "LGALS7_detection_by_dataset.csv"), row.names = FALSE)
  plot_df <- focus %>% filter(n_cells >= 30)
  if (nrow(plot_df) > 0) {
    plot_df$cell_type <- factor_celltype(plot_df$cell_type)
    p <- ggplot(plot_df, aes(cell_type, pct_detected, fill = cell_type)) +
      geom_col() +
      facet_grid(dataset ~ tissue, scales = "free_x") +
      scale_fill_manual(values = colors_for(levels(plot_df$cell_type))) +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(
        x = NULL,
        y = paste0("% cells with ", FOCUS_GENE, " detected"),
        title = "Galectin-7 detection across cohorts",
        subtitle = "Percent detected is comparable across datasets; mean expression is not"
      )
    save_plot(p, file.path(out, "LGALS7_detection_across_datasets"), width = 12, height = 8)
  }
  invisible(focus)
}

main <- function() {
  message("Project: ", PROJECT_DIR)
  message("Datasets: ", paste(DATASETS, collapse = ", "))
  message(
    "Local IMC files are not used. The 25-protein panel has no galectin-7 measurement."
  )
  errors <- list()
  for (dataset in DATASETS) {
    result <- tryCatch(
      {
        analyze_dataset(dataset)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
    if (!is.null(result)) {
      errors[[dataset]] <- result
      err_dir <- file.path(OUT_DIR, dataset)
      dir.create(err_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines(result, file.path(err_dir, "ERROR.txt"))
      message("FAILED ", dataset, ": ", result)
    }
  }
  combine_datasets()
  writeLines(
    capture.output(sessionInfo()),
    file.path(OUT_DIR, "sessionInfo.txt")
  )
  if (length(errors) > 0) {
    stop(
      "One or more datasets failed:\n",
      paste(names(errors), errors, sep = ": ", collapse = "\n"),
      call. = FALSE
    )
  }
  message("Done. Results are in ", OUT_DIR)
}

if (sys.nframe() == 0L) {
  main()
}
