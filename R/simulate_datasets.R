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
#' processing times
simCompDataset <- function(dataset, root, seed = 1,
                           bp = BiocParallel::SerialParam()) {

    timings <- matrix(nrow = 7, ncol = 2)
    rownames(timings) <- c("Splat", "SplatDrop", "Simple", "Lun", "Lun2",
                           "Lun2ZINB", "scDD")
    colnames(timings) <- c("Estimation", "Simulation")

    set.seed(seed)
    message("#### STARTING ", dataset["Dataset"], " ####")
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
    tt <- system.time(params <- splatter::splatEstimate(counts))
    timings["Splat", "Estimation"] <- tt[3]
    timings["SplatDrop", "Estimation"] <- tt[3]

    timings["Splat", "Simulation"] <- system.time(
    sims$Splat <- splatter::splatSimulate(params, dropout.present = FALSE,
                                          verbose = FALSE, seed = seed)
    )[3]
    timings["SplatDrop", "Simulation"] <- system.time(
    sims$SplatDrop <- splatter::splatSimulate(params, dropout.present = TRUE,
                                              verbose = FALSE, seed = seed)
    )[3]

    message("***ADDING SIMPLE***")
    timings["Simple", "Estimation"] <- system.time(
    params <- splatter::simpleEstimate(counts)
    )[3]
    timings["Simple", "Simulation"] <- system.time(
    sims$Simple <- splatter::simpleSimulate(params, verbose = FALSE,
                                            seed = seed)
    )[3]

    message("***ADDING LUN***")
    timings["Lun", "Estimation"] <- system.time(
    params <- splatter::lunEstimate(counts)
    )[3]
    timings["Lun", "Simulation"] <- system.time(
    sims$Lun <- splatter::lunSimulate(params, verbose = FALSE, seed = seed)
    )[3]

    message("***ADDING LUN2***")
    tt <- system.time(
    params <- splatter::lun2Estimate(counts,
                                     plate = sample(1:2, ncol(counts),
                                                    replace = TRUE),
                                     min.size = 20, BPPARAM = bp,
                                     verbose = FALSE)
    )
    timings["Lun2", "Estimation"] <- tt[3]
    timings["Lun2ZINB", "Estimation"] <- tt[3]
    timings["Lun2", "Simulation"] <- system.time(
    sims$Lun2 <- splatter::lun2Simulate(params, nGenes = nrow(counts),
                                        verbose = FALSE, seed = seed)
    )[3]
    timings["Lun2ZINB", "Simulation"] <- system.time(
    sims$Lun2ZINB <- splatter::lun2Simulate(params, nGenes = nrow(counts),
                                            zinb = TRUE, verbose = FALSE,
                                            seed = seed)
    )[3]

    message("***ADDING SCDD***")
    timings["scDD", "Estimation"] <- system.time(
    params <- splatter::scDDEstimate(counts,
                                     conditions = sample(1:2, ncol(counts),
                                                         replace = TRUE))
    )[3]

    half.nGenes <- round(nrow(counts) / 2)
    timings["scDD", "Simulation"] <- system.time(
    sim.scDD <- splatter::scDDSimulate(params, nCells = round(ncol(counts) / 2),
                                       nDE = 0, nDP = 0, nDM = 0, nDB = 0,
                                       nEE = half.nGenes,
                                       nEP = nrow(counts) - half.nGenes,
                                       seed = seed, verbose = FALSE,
                                       BPPARAM = bp)
    )[3]
    sims$scDD <- sim.scDD

    message("***TESTING GENE GoF***")
    sims <- bplapply(sims, testGenesGoF, BPPARAM = bp)

    message("***COMPARING DATASETS***")
    cols <- scales::hue_pal()(length(sims))
    comp <- splatter::compareSCESets(sims, point.size = 0.3, colours = cols)
    diff <- splatter::diffSCESets(sims, ref = "Real", point.size = 0.3,
                                  colours = cols[-1])

    message("#### DONE ", sum(timings), "####")

    return(list(Comp = comp, Diff = diff, Timings = timings))
}
