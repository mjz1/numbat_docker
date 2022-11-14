#!/usr/bin/env Rscript

####
# The default docker doesn't have this function. It's the only external dependency so we just copy the seurat function here
Read10X <- function(
  data.dir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
) {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, 'barcodes.tsv')
    gene.loc <- file.path(run, 'genes.tsv')
    features.loc <- file.path(run, 'features.tsv.gz')
    matrix.loc <- file.path(run, 'matrix.mtx')
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc) ) {
      stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    } else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names,
        FUN = ExtractField,
        field = 1,
        delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(
      file = ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc),
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning(
        'Some features names are NA. Replacing NA names with ID from the opposite column requested',
        call. = FALSE,
        immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column,
                    " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                    " Try setting the gene.column argument to a value <= to ", fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added
    # and we will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) { # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(
        X = lvls,
        FUN = function(l) {
          return(data[data_types == l, , drop = FALSE])
        }
      )
      names(x = data) <- lvls
    } else{
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    # list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
    list_of_data[[j]] <- as(list_of_data[[j]], "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}
#####

# Script starts here
library(argparse)

# Parse command line arguments
parser <- ArgumentParser()

parser$add_argument("--patient", help="Patient identifier")
parser$add_argument("--samples", help="Patient sample names")
parser$add_argument("--mtxdirs", help="comma seperated list of matrix folders")
parser$add_argument("--outdir", help="output directory")
parser$add_argument("--genome_ver", help="Genome version (hg19, hg38) [default %(default)s]", default = "hg38")
parser$add_argument("--cores", help="number of cores to use [default %(default)s]", default=4, type = "integer")
parser$add_argument("--trans", help="Numbat HMM transmission probability [default %(default)s]", default=1e-5, type="double")
parser$add_argument("--gamma", help="Numbat overdispersion parameter in allele counts [default %(default)s]", default=20, type = "integer")
parser$add_argument("--min_cells", help="Numbat minimum number of cells for which an pseudobulk HMM will be run [default %(default)s]", default=20, type="integer")
parser$add_argument("--min_LLR", help="Numbat minimum log-likelihood ratio threshold to filter CNVs by. [default %(default)s]", default=50, type="integer")
parser$add_argument("--multi_allelic", help="Flag to run numbat in multi-allelic mode [default %(default)s]", action="store_true", default=TRUE)
parser$add_argument("--high_purity", help="Flag to detect and exclude regions of clonal deletions/LOH before running Numbat. Recommended for cell line data or high-purity tumors", action="store_true", default=FALSE)

args <- parser$parse_args()

patient = args$patient
samples = args$samples
mtxdirs = args$mtxdirs
outdir = args$outdir
genome_ver = args$genome_ver
ncores = args$cores
trans = args$trans
gamma = args$gamma
min_cells = args$min_cells
min_LLR = args$min_LLR
multi_allelic = args$multi_allelic
high_purity = args$high_purity

setwd(outdir)

# Load libraries
library(numbat)
# library(Seurat)

# For multi-sample patients parse and load
samps <- strsplit(samples, ",")[[1]]
mtx_dirs <- strsplit(mtxdirs, ",")[[1]]

allele_files <- file.path(outdir, paste0(samps, "_allele_counts.tsv.gz"))

if (length(samps) > 1) {
    message("Multiple samples detected. Reading and merging data...")    
}

df_allele <- c()
count_mat <- c()
for (i in 1:length(samps)) {
    # Load count matrix and relabel cells by their sample
    count_mat_ <- Read10X(mtx_dirs[i])

    # Check for multiomic and only return gene expression
    if (length(count_mat_) > 1 & !is.null(names(count_mat_))) {
        print("Multiomic data detected. Keeping only gene expression...")
        count_mat_ <- count_mat_[['Gene Expression']]
    }

    # Paste sample names into cell_ids
    colnames(count_mat_) <- paste(samps[i], "_", colnames(count_mat_), sep = "")
    count_mat <- cbind(count_mat, count_mat_)

    # Read and add cols to the allele df
    df_allele_ <- read.table(file = allele_files[i], header = T, sep = "\t")
    df_allele_$cell <- paste(samps[i], "_", df_allele_$cell, sep = "")
    df_allele_$sample <- samps[i]
    df_allele_$group <- patient
    df_allele <- rbind(df_allele, df_allele_)
}


# Save input files
save(x = count_mat, file = "count_mat.rda")
save(x = df_allele, file = "df_allele.rda")

# Set multi_allelic flag
if (multi_allelic) {
    multi_allelic=TRUE
} else {
    multi_allelic=FALSE
}

# Detect clonal loh for high purity flagged samples
segs_loh = NULL
if (high_purity) {
    message("Running with high purity flag on!")
    bulk = get_bulk(
        count_mat = count_mat,
        lambdas_ref = ref_hca,
        df_allele = df_allele,
        gtf = gtf_hg38
    )

    segs_loh = detect_clonal_loh(bulk, t = 1e-4)

    # Save a file to keep track of high purity mode runs
    fileConn<-file("high_purity")
    writeLines(c(""), fileConn)
    close(fileConn)
}

# Run numbat
out = run_numbat(
    count_mat = count_mat, # gene x cell integer UMI count matrix 
    lambdas_ref = ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
    df_allele = df_allele, # allele dataframe generated by pileup_and_phase script
    genome = genome_ver,
    min_cells = min_cells,
    t = trans,
    max_iter = 2,
    min_LLR = min_LLR,
    init_k = 3,
    ncores = ncores,
    plot = TRUE,
    out_dir = outdir,
    multi_allelic = multi_allelic,
    segs_loh = segs_loh
)

save(x = out, file = "numbat_out.rda")
