#' Test gene goodness-of-fit
#'
#' Test the goodness of fit for each gene with regards to the negative-binomial,
#' log-normal and normal distributions.
#'
#' @param sce SCESet object
#'
#' @details For each gene the parameters for the distributions are estimated
#' using \code{\link[fitdistrplus]{fitdist}} then the counts for the gene are
#' tested against a theoretical distribution with those parameters using a
#' Chi-squared goodness-of-fit test.
#'
#' @return SCESet object with additional fData columns
testGenesGoF <- function(sce) {
    counts <- scater::counts(sce)

    stats <- apply(counts, 1, function(gene.counts) {
        fit <- try(fitdistrplus::fitdist(gene.counts, "nbinom"), silent = TRUE)

        if (class(fit) != "try-error") {
            mu <- fit$estimate["mu"]
            size <- fit$estimate["size"]

            gene.hist <- hist(gene.counts, breaks = 10, right = FALSE,
                              plot = FALSE)

            breaks.cdf <- pnbinom(gene.hist$breaks, mu = mu, size = size)
            null.probs <- zoo::rollapply(breaks.cdf, 2,
                                         function(x) {x[2] - x[1]})

            test <- chisq.test(gene.hist$counts, p = null.probs,
                               rescale.p = TRUE)

            out.nb <- c(test$statistic, test$p.value)
        } else {
            out.nb <- c(NA, NA)
        }

        fit <- try(fitdistrplus::fitdist(gene.counts + 1, "lnorm"),
                   silent = TRUE)

        if (class(fit) != "try-error") {
            loc <- fit$estimate["meanlog"]
            scale <- fit$estimate["sdlog"]

            gene.hist <- hist(gene.counts + 1, breaks = 10, right = FALSE,
                              plot = FALSE)

            breaks.cdf <- plnorm(gene.hist$breaks, meanlog = loc, sdlog = scale)
            null.probs <- zoo::rollapply(breaks.cdf, 2,
                                         function(x) {x[2] - x[1]})

            test <- chisq.test(gene.hist$counts, p = null.probs,
                               rescale.p = TRUE)

            out.ln <- c(test$statistic, test$p.value)
        } else {
            out.ln <- c(NA, NA)
        }

        fit <- try(fitdistrplus::fitdist(gene.counts, "norm"), silent = TRUE)

        if (class(fit) != "try-error") {
            mean <- fit$estimate["mean"]
            sd <- fit$estimate["sd"]

            gene.hist <- hist(gene.counts, breaks = 10, right = FALSE,
                              plot = FALSE)

            breaks.cdf <- plnorm(gene.hist$breaks, mean = mean, sd = sd)
            null.probs <- zoo::rollapply(breaks.cdf, 2,
                                         function(x) {x[2] - x[1]})

            test <- chisq.test(gene.hist$counts, p = null.probs,
                               rescale.p = TRUE)

            out.norm <- c(test$statistic, test$p.value)
        } else {
            out.norm <- c(NA, NA)
        }

        out <- c(out.nb, out.ln, out.norm)
        names(out) <- c("NBChi", "NBPVal", "LNChi", "LNPVal", "NormChi",
                        "NormPVal")

        return(out)
    })

    stats <- data.frame(t(stats))

    Biobase::fData(sce) <- cbind(Biobase::fData(sce), stats)

    return(sce)
}
