#' Load Tung
#'
#' Load various versions of the Tung data and associated annotation
#'
#' @param count The type of counts to load. Either "reads" for full length reads
#'        or "molecules" for UMI counts.
#' @param version The version of the data to load. Choice of "raw" (all genes
#'        and cells), "nozeros" (all zero genes removed), "filter" (bad cells
#'        and genes removed) or "final" (UMI expression data after final
#'        processing).
#' @param individual Individual to select. Options are "all" (default),
#'        "NA19098", "NA19101" or "NA19239".
#' @param bad.cell Vector of "bad" cell names.
#' @param remove.bad Whether to remove "bad" cells from "raw" or "nozeros" data.
#' @param remove.ercc Whether to remove features corresponding to the ERCC
#'        spike-ins.
#' @param remove.zeros Whether to remove all-zero genes.
#'
#' @return SCESet containing the chosen dataset and annotation
#' @examples
#' # Load raw read counts
#' tung <- loadTung(count = "reads", version = "raw")
#' # Load filtered UMI counts
#' tung <- loadTung(count = "molecules", version = "filter")
loadTung <- function(count = c("reads", "molecules"),
                     version = c("raw", "nozeros", "filter", "final"),
                     individual = c("all", "NA19098", "NA19101", "NA19239"),
                     bad.cells = NULL, remove.bad = FALSE, remove.ercc = FALSE,
                     remove.zeros = FALSE) {

    count <- match.arg(count)
    version <- match.arg(version)
    individual <- match.arg(individual)

    # Choose file to load
    dir <- "/group/bioi1/shared/public_data/GiladSingleCellHapMap/data/"

    count.cols <- readr::cols(.default = readr::col_integer(),
                              X1 = readr::col_character())
    switch(version,
           raw = {
               count.file <- paste0(count, "-raw-single-per-sample.txt")
               count.cols <- readr::cols(
                   .default = readr::col_integer(),
                   individual = readr::col_character(),
                   replicate = readr::col_character(),
                   well = readr::col_character()
               )
           },
           nozeros = {
               count.file <- paste0(count, ".txt")
           },
           filter = {
               count.file <- paste0(count, "-filter.txt")
           },
           final = {
               count.file <- paste0(count, "-final.txt")
           }
    )

    if ((count == "reads") & (version == "final")) {
        stop("There is no 'final' version for reads. Please try again.")
    }

    if (remove.bad && is.null(bad.cells)) {
        warning("No bad cells supplied. No cells will be removed.")
        remove.bad = FALSE
    }

    # Read data
    message("Reading '", count.file, "'...")
    counts.raw <- readr::read_tsv(file.path(dir, count.file),
                                  col_types = count.cols)

    if (version == "raw") {
        # Raw data is transposed and has extra information
        samples <- paste(counts.raw$individual, counts.raw$replicate,
                         counts.raw$well, sep = ".")
        counts <- t(counts.raw[, -c(1:3)])
        colnames(counts) <- samples
        genes <- rownames(counts)
    } else {
        genes <- counts.raw[[1]]
        counts <- as.matrix(counts.raw[, -1])
        rownames(counts) <- genes
        samples <- colnames(counts)
    }

    # Get ERCC gene names
    genes.ercc <- grep(pattern = "ERCC-", genes, value = TRUE)

    message("Getting phenotype data...")
    annot <- readr::read_tsv(file.path(dir, "annotation.txt"),
                             col_types = readr::cols(
                                 .default = readr::col_character()))
    qc <- readr::read_tsv(file.path(dir, "qc-ipsc.txt"),
                          col_types = readr::cols(
                              .default = readr::col_character(),
                              cell_number = readr::col_integer(),
                              concentration = readr::col_double(),
                              tra1.60 = readr::col_integer()
                          ))
    qc <- dplyr::mutate(qc, sample = paste(individual, replicate, well,
                                           sep = "."))
    qc <- qc[, -c(1:3)]
    pheno.data <- annot[, c(5, 1:4)]
    pheno.data <- dplyr::left_join(pheno.data, qc, by = c(sample_id = "sample"))
    if (!is.null(bad.cells)) {
        pheno.data <- dplyr::mutate(pheno.data,
                                    bad_cell = sample_id %in% bad.cells)
    }
    # SCESet needs dataframe with correct sample names in order
    pheno.data <- as.data.frame(pheno.data)
    rownames(pheno.data) <- pheno.data$sample_id
    pheno.data <- pheno.data[samples, ]

    message("Getting feature data...")
    gene.info.fc <- readr::read_tsv(file.path(dir, "count-matrix.txt"),
                                    comment = "#",
                                    col_types = readr::cols(
                                        .default = readr::col_integer(),
                                        Geneid = readr::col_character(),
                                        Chr = readr::col_character(),
                                        Strand = readr::col_character()))[, 1:6]
    gene.info <- readr::read_tsv(file.path(dir, "gene-info.txt"),
                                 col_types = readr::cols(
                                     .default = readr::col_character(),
                                     transcript_count = readr::col_integer()))
    gene.info <- gene.info[, -2]
    cc.genes <- readr::read_tsv(file.path(dir, "cellcyclegenes.txt"),
                                col_types = readr::cols(
                                    .default = readr::col_character()))
    pp.genes <- readr::read_tsv(file.path(dir, "pluripotency-genes.txt"),
                                col_types = readr::cols(
                                    .default = readr::col_character()))
    feat.data <- dplyr::filter(gene.info.fc, Geneid %in% genes)
    feat.data <- dplyr::left_join(feat.data, gene.info,
                                  by = c(Geneid = "ensembl_gene_id"))
    # Add stages of cell cycle as columns
    feat.data <- dplyr::mutate(feat.data, "cc_G1/S" = Geneid %in% cc.genes$`G1/S`)
    feat.data <- dplyr::mutate(feat.data, "cc_S" = Geneid %in% cc.genes$S)
    feat.data <- dplyr::mutate(feat.data, "cc_G2/M" = Geneid %in% cc.genes$`G2/M`)
    feat.data <- dplyr::mutate(feat.data, "cc_M" = Geneid %in% cc.genes$M)
    feat.data <- dplyr::mutate(feat.data, "cc_M/G1" = Geneid %in% cc.genes$`M/G1`)
    feat.data <- dplyr::mutate(feat.data, cc_any = `cc_G1/S` | `cc_S` |
                                                   `cc_G2/M` | `cc_M` |
                                                   `cc_M/G1`)
    feat.data <- dplyr::mutate(feat.data, pluripotent = Geneid %in% pp.genes$To)
    # SCESet need dataframe with correct gene names in order
    feat.data <- as.data.frame(feat.data)
    rownames(feat.data) <- feat.data$Geneid
    feat.data <- feat.data[genes, ]

    # Produce SCESet object
    message("Building SCESet object...")
    if (version == "final") {
        sce <- scater::newSCESet(exprsData = counts)
    } else {
        sce <- scater::newSCESet(countData = counts)
    }
    message("Adding annotation...")
    Biobase::pData(sce) <- pheno.data
    Biobase::fData(sce) <- feat.data
    if (remove.bad) {
        message("Removing bad cells...")
        sce <- sce[, !(Biobase::pData(sce)$sample_id %in% bad.cells)]
    }
    if (individual != "all") {
        message("Selecting individual...")
        sce <- sce[, Biobase::pData(sce)$individual == individual]
    }
    if (remove.ercc) {
        message("Removing ERCC features...")
        sce <- sce[!(Biobase::fData(sce)$Geneid %in% genes.ercc), ]
    }
    if (remove.zeros) {
        message("Removing all-zero features...")
        sce <- sce[!(rowSums(scater::counts(sce)) == 0), ]
    }
    if (version != "final") {
        message("Calculating TPM...")
        lengths <- Biobase::fData(sce)$Length
        scater::tpm(sce) <- scater::calculateTPM(sce,
                                                 effective_length = lengths)
    }
    message("Calculating QC metrics...")
    sce <- scater::calculateQCMetrics(sce, feature_controls = genes.ercc)

    message("Done!")
    return(sce)
}
