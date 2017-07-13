#' Test gene goodness-of-fit
#'
#' Test the goodness of fit for each gene with regards to the negative-binomial,
#' log-normal and normal distributions.
#'
#' @param sce SCESet object
#'
#' @details For each gene the parameters for the distributions are estimated
#' using \code{\link[fitdistrplus]{fitdist}} then the goodness of fit is tested
#' using the Chi-squared statistic from \code{\link[fitdistrplus]{gofstat}}.
#'
#' @return SCESet object with additional fData columns
testGenesGoF <-  function(sce) {
    counts <- scater::counts(sce)

    dists <- c(NB = "nbinom", LN = "lnorm", Norm = "norm", Poi = "pois")

    stats <- apply(counts, 1, function(gene.counts) {

        outs <- list()

        for (dist in dists) {
            dummy <- capture.output(
            fit <- try(fitdistrplus::fitdist(gene.counts, dist), silent = TRUE)
            )

            out <- c(NA, NA)
            if (class(fit) != "try-error") {
                gof <- try(fitdistrplus::gofstat(fit), silent = TRUE)

                if (class(gof) != "try-error") {
                    out <- c(gof$chisq, gof$chisqpvalue)
                }
            }

            outs <- c(outs, out)
        }

        outs <- do.call(c, outs)
    })

    stats <- data.frame(t(stats))
    colnames(stats) <- paste0(rep(names(dists), each = 2), c("Chi", "PVal"))

    Biobase::fData(sce) <- cbind(Biobase::fData(sce), stats)

    return(sce)
}
