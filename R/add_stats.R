addGeneStats <- function(sce, value = c("counts", "cpm", "tpm", "fpkm"),
                         log = FALSE, offset = 1, no.zeros = FALSE) {

    value <- match.arg(value)

    switch(value,
           counts = {
               values = scater::counts(sce)
           },
           cpm = {
               values = scater::cpm(sce)
           },
           tpm = {
               values = scater::tpm(sce)
           },
           fpkm = {
               values = scater::fpkm(sce)
           }
    )

    suffix = value

    if (no.zeros) {
        values[values == 0] <- NA
        suffix = paste0(suffix, "_no0")
    }

    if (log) {
        values = log2(values + offset)
        suffix = paste0("log_", suffix)
    }

    mean.str <- paste0("mean_", suffix)
    var.str  <- paste0("var_",  suffix)
    cv.str   <- paste0("cv_",   suffix)
    med.str  <- paste0("med_",  suffix)
    mad.str  <- paste0("mad_",  suffix)

    Biobase::fData(sce)[, mean.str] <- rowMeans(values, na.rm = TRUE)
    Biobase::fData(sce)[, var.str]  <- matrixStats::rowVars(values,
                                                             na.rm = TRUE)
    Biobase::fData(sce)[, cv.str]   <- sqrt(Biobase::fData(sce)[, var.str]) /
        Biobase::fData(sce)[, mean.str]
    Biobase::fData(sce)[, med.str]  <- matrixStats::rowMedians(values,
                                                                na.rm = TRUE)
    Biobase::fData(sce)[, mad.str]  <- matrixStats::rowMads(values,
                                                             na.rm = TRUE)

    return(sce)
}

addIsZeroOutlier <- function(sce) {
    med <- Biobase::fData(sce)$med_log_cpm
    mad <- Biobase::fData(sce)$mad_log_cpm
    Biobase::fData(sce)$is_zero_outlier <- (med / mad) > 3

    return(sce)
}
