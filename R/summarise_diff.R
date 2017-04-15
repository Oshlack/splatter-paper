#' Summarise diffSCESets
#'
#' Summarise the results of Splatter's diffSCESets function. The various
#' properties are sorted, differences calculated, the Median Absolute Deviation
#' taken as the statistic and the ranks calculated.
#'
#' @param diff Output of Spltter's diffSCESets function
#'
#' @return Summary table with MADs, ranks and both combined in long format
summariseDiff <- function(diff) {

    dataset.names <- unique(diff$PhenoData$Dataset)

    # Calculate gene level MADs
    fData.mads <- sapply(dataset.names, function(dataset) {
        df <- diff$FeatureData[diff$FeatureData$Dataset == dataset, ]
        mean <- median(abs(df$RankDiffMeanLogCPM))
        var <- median(abs(df$RankDiffVarLogCPM))
        zeros <- median(abs(df$RankDiffZeros))
        mean.var <- median(abs(df$MeanRankVarDiff))
        mean.zeros <- median(abs(df$MeanRankZerosDiff))
        return(c(Mean = mean, Variance = var, ZerosGene = zeros,
                 MeanVar = mean.var, MeanZeros = mean.zeros))
    })

    # Calculate cell level MADs
    pData.mads <- sapply(dataset.names, function(dataset) {
        df <- diff$PhenoData[diff$PhenoData$Dataset == dataset, ]
        lib.size <- median(abs(df$RankDiffLibSize))
        zeros <- median(abs(df$RankDiffZeros))
        return(c(LibSize = lib.size, ZerosCell = zeros))
    })

    # Combine them
    mads <- data.frame(Dataset = dataset.names, t(fData.mads), t(pData.mads))

    # Calculate the ranks
    fData.ranks <- matrixStats::rowRanks(fData.mads)
    pData.ranks <- matrixStats::rowRanks(pData.mads)

    ranks <- data.frame(Dataset = dataset.names, t(fData.ranks), t(pData.ranks))
    colnames(ranks) <- paste0(colnames(mads), "Rank")

    # Convert to long format (for easier ggplot plotting)
    mads.long <- stats::reshape(mads, varying = 2:8, direction = "long",
                                idvar = "Dataset", timevar = "Statistic",
                                times = colnames(mads)[2:8], v.names = "MAD")
    ranks.long <- stats::reshape(ranks, varying = 2:8, direction = "long",
                                 idvar = "Dataset", timevar = "Statistic",
                                 times = colnames(ranks)[2:8], v.names = "Rank")

    # Combine the long data
    long <- data.frame(mads.long, Rank = ranks.long$Rank)
    row.names(long) <- NULL

    # Combine MADs, ranks and long data
    summary <- list(MADs = mads, Ranks = ranks, Long = long)

    return(summary)
}
