#' Simulate and compare dataset
#'
#' Takes a vector describing a dataset, loads the datasets, estimates parameters
#' and simulates data using various models and retuns a comparison.
#'
#' @param dataset Vector describing a real dataset
#' @param root Root path where data is stored
#' @param seed Random seed
#' @param bp BiocParallel BPPARAM object
#'
#' @return List with the results of Splatter's comparison functions and the
#' processin time
simCompDataset <- function(dataset, root, seed = 1,
                           bp = BiocParallel::SerialParam()) {
    tt <- system.time({ # Start timed section
    set.seed(seed)
    message("#### STARTING ", dataset["Dataset"], "####")
    message("***LOADING COUNTS***")
    counts <- loadDataset(dataset, root)
    counts <- counts[rowSums(counts) > 0, ]
    na.rows <- which(rowSums(is.na(counts)) > 0)
    if (length(na.rows) > 0) {
        counts <- counts[-na.rows, ]
    }
    counts <- counts[, sample(1:ncol(counts), 200)]

    sims <- list(Real = scater::newSCESet(countData = counts))

    message("***ADDING SPLAT***")
    params <- splatter::splatEstimate(counts)
    sims$Splat <- splatter::splatSimulate(params, dropout.present = FALSE,
                                          verbose = FALSE, seed = seed)
    sims$SplatDrop <- splatter::splatSimulate(params, dropout.present = TRUE,
                                              verbose = FALSE, seed = seed)

    message("***ADDING SIMPLE***")
    params <- splatter::simpleEstimate(counts)
    sims$Simple <- splatter::simpleSimulate(params, verbose = FALSE,
                                            seed = seed)

    message("***ADDING LUN***")
    params <- splatter::lunEstimate(counts)
    sims$Lun <- splatter::lunSimulate(params, verbose = FALSE, seed = seed)

    message("***ADDING LUN2***")
    params <- splatter::lun2Estimate(counts,
                                     plate = sample(1:2, ncol(counts),
                                                    replace = TRUE),
                                     min.size = 20, BPPARAM = bp,
                                     verbose = FALSE)
    sims$Lun2 <- splatter::lun2Simulate(params, nGenes = nrow(counts),
                                        verbose = FALSE, seed = seed)
    sims$Lun2ZINB <- splatter::lun2Simulate(params, nGenes = nrow(counts),
                                            zinb = TRUE, verbose = FALSE,
                                            seed = seed)

    message("***ADDING SCDD***")
    params <- splatter::scDDEstimate(counts,
                                     conditions = sample(1:2, ncol(counts),
                                                         replace = TRUE))

    half.nGenes <- round(nrow(counts) / 2)
    sim.scDD <- splatter::scDDSimulate(params, nCells = round(ncol(counts) / 2),
                                       nDE = 0, nDP = 0, nDM = 0, nDB = 0,
                                       nEE = half.nGenes,
                                       nEP = nrow(counts) - half.nGenes,
                                       seed = seed, verbose = FALSE,
                                       BPPARAM = bp)
    sims$scDD <- sim.scDD

    cols <- scales::hue_pal()(length(sims))

    message("***COMPARING DATASETS***")
    comp <- splatter::compareSCESets(sims, point.size = 0.3, colours = cols)
    diff <- splatter::diffSCESets(sims, ref = "Real", point.size = 0.3,
                                  colours = cols[-1])
    })[3] # End timed section

    message("#### DONE ", tt, "####")

    return(list(Comp = comp, Diff = diff, Time = tt))
}
