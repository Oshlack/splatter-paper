#' Simulate dataset
#'
#' Takes a counts matrix, estimates parameters and simulates data using various
#' models.
#'
#' @param counts Counts matrix to estimate parameters from
#' @param sims Vector of simulations to compare
#' @param seed Random seed
#' @param bp BiocParallel BPPARAM object
#'
#' @return List with the simulations and the processing times
simData <- function(counts, models = c("Splat", "SplatDrop", "Simple", "Lun",
                                       "Lun2", "Lun2ZINB", "scDD", "BASiCS"),
                    seed = 1, verbose = FALSE,
                    bp = BiocParallel::SerialParam()){

    # Total processing time
    ttime <- system.time({

    timings <- matrix(nrow = length(models), ncol = 2)
    rownames(timings) <- models
    colnames(timings) <- c("Estimation", "Simulation")

    set.seed(seed)

    sims <- list(Real = scater::newSCESet(countData = counts))

    if ("Splat" %in% models || "SplatDrop" %in% models) {
        message("***ADDING SPLAT***")
        tt <- system.time(params <- splatter::splatEstimate(counts))
    }

    if ("Splat" %in% models) {
        timings["Splat", "Estimation"] <- tt[3]
        timings["Splat", "Simulation"] <- system.time(
            sims$Splat <- splatter::splatSimulate(params,
                                                  dropout.present = FALSE,
                                                  verbose = verbose,
                                                  seed = seed)
        )[3]
    }

    if ("SplatDrop" %in% models) {
        timings["SplatDrop", "Estimation"] <- tt[3]
        timings["SplatDrop", "Simulation"] <- system.time(
            sims$SplatDrop <- splatter::splatSimulate(params,
                                                      dropout.present = TRUE,
                                                      verbose = verbose,
                                                      seed = seed)
        )[3]
    }

    if ("Simple" %in% models) {
        message("***ADDING SIMPLE***")
        timings["Simple", "Estimation"] <- system.time(
            params <- splatter::simpleEstimate(counts)
        )[3]
        timings["Simple", "Simulation"] <- system.time(
            sims$Simple <- splatter::simpleSimulate(params, verbose = verbose,
                                                    seed = seed)
        )[3]
    }


    if ("Lun" %in% models) {
        message("***ADDING LUN***")
        timings["Lun", "Estimation"] <- system.time(
            params <- splatter::lunEstimate(counts)
        )[3]
        timings["Lun", "Simulation"] <- system.time(
            sims$Lun <- splatter::lunSimulate(params, verbose = verbose,
                                              seed = seed)
        )[3]
    }

    if ("Lun2" %in% models || "Lun2ZINB" %in% models) {
        message("***ADDING LUN2***")
        tt <- system.time(
            params <- splatter::lun2Estimate(counts,
                                             plate = sample(1:2, ncol(counts),
                                                            replace = TRUE),
                                             min.size = 20, BPPARAM = bp,
                                             verbose = verbose)
        )
    }

    if ("Lun2" %in% models) {
        timings["Lun2", "Estimation"] <- tt[3]
        timings["Lun2", "Simulation"] <- system.time(
            sims$Lun2 <- splatter::lun2Simulate(params, nGenes = nrow(counts),
                                                verbose = verbose, seed = seed)
        )[3]
    }

    if ("Lun2ZINB" %in% models) {
        timings["Lun2ZINB", "Estimation"] <- tt[3]
        timings["Lun2ZINB", "Simulation"] <- system.time(
            sims$Lun2ZINB <- splatter::lun2Simulate(params,
                                                    nGenes = nrow(counts),
                                                    zinb = TRUE,
                                                    verbose = verbose,
                                                    seed = seed)
        )[3]
    }

    if ("scDD" %in% models) {
        message("***ADDING SCDD***")
        timings["scDD", "Estimation"] <- system.time(
            params <- splatter::scDDEstimate(counts,
                                             conditions = sample(1:2, ncol(counts),
                                                                 replace = TRUE),
                                             BPPARAM = bp)
        )[3]

        timings["scDD", "Simulation"] <- system.time(
            sims$scDD <- splatter::scDDSimulate(params,
                                                nCells = round(ncol(counts) / 2),
                                                seed = seed, verbose = verbose,
                                                BPPARAM = bp)
        )[3]
    }

    if ("BASiCS" %in% models) {
        message("***ADDING BASiCS***")
        timings["BASiCS", "Estimation"] <- system.time(
            params <- splatter::BASiCSEstimate(counts,
                                               batch = sample(1:2, ncol(counts),
                                                              replace = TRUE),
                                               verbose = verbose,
                                               progress = verbose)
        )[3]

        timings["BASiCS", "Simulation"] <- system.time(
            sims$BASiCS <- splatter::BASiCSSimulate(params, seed = seed,
                                                    verbose = verbose)
        )[3]
    }

    })

    message("#### DONE ", ttime[3] / 60 / 60, " hours ####")

    return(list(Sims = sims, Timings = timings))
}

#' @param compare Whether to produce a comparison of simulations
#' @param test.gof Whether to test gene goodness-of-fit
#if (test.gof) {
#    message("***TESTING GENE GoF***")
#    sims <- BiocParallel::bplapply(sims, testGenesGoF, BPPARAM = bp)
#}
#
#if (compare) {
#    message("***COMPARING DATASETS***")
#    cols <- scales::hue_pal()(length(sims))
#    comp <- splatter::compareSCESets(sims, point.size = 0.3, colours = cols)
#    diff <- splatter::diffSCESets(sims, ref = "Real", point.size = 0.3,
#                                  colours = cols[-1])
#} else {
#    comp = NA
#    diff = NA
#}
