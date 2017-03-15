loadDataset <- function(dataset, root) {

    file <- file.path(root, dataset["Path"], dataset["CountsFile"])

    if (dataset["BpipePipeline"] == "Yes") {
        counts <- readr::read_tsv(file, comment = "#",
                                  col_types = readr::cols(
                                      .default = readr::col_integer(),
                                      Geneid = readr::col_character(),
                                      Chr = readr::col_character(),
                                      Start = readr::col_character(),
                                      End = readr::col_character(),
                                      Strand = readr::col_character()
                                  )
        )

        counts <- dplyr::select(counts, -Geneid, -Chr, -Start, -End, -Strand)
        counts <- as.matrix(counts)

    } else if (dataset["Conquer"] == "Yes") {
        counts <- readr::read_tsv(file,
                                  col_types = readr::cols(
                                      .default = readr::col_double(),
                                      Gene = readr::col_character()
                                  )
        )

        counts <- dplyr::select(counts, -Gene)
        counts <- as.matrix(counts)

    } else if (dataset["Dataset"] == "Grun") {
        counts <- readr::read_tsv(file,
                                  col_types = readr::cols(
                                      .default = readr::col_double(),
                                      GENENAME = readr::col_character()
                                      )
                                  )

        counts <- dplyr::select(counts, -GENENAME)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Klein") {
        counts <- readr::read_csv(file,
                                  col_types = readr::cols(
                                      .default = readr::col_integer(),
                                      Gene = readr::col_character()
                                      ),
                                  col_names = c("Gene", paste0("S", 1:239)),
                                  skip = 1)

        counts <- dplyr::select(counts, -Gene)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Tung") {
        counts <- readr::read_tsv(file,
                                  col_types = readr::cols(
                                      .default = readr::col_integer(),
                                      X1 = readr::col_character()
                                  ))

        counts <- dplyr::select(counts, -X1)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Zeisel") {
        counts <- read_tsv(file,
                           col_types = readr::cols(
                               .default = readr::col_character()
                               ),
                           col_names = FALSE,
                           skip = 11)

        counts <- dplyr::select(counts, -X1)
        counts <- as.matrix(counts)
        counts <- apply(counts, 1, as.numeric)
    } else if (dataset["Dataset"] == "Zieg") {
        counts <- readr::read_tsv(file,
                                  col_types = readr::cols(
                                      .default = readr::col_integer(),
                                      Gene = readr::col_character()
                                  ),
                                  col_names = c("Gene", paste0("S", 1:482)),
                                  skip = 1)
        counts <- dplyr::select(counts, -Gene)
        counts <- as.matrix(counts)
    } else if (dataset["Dataset"] == "Velten") {
        counts <- readr::read_csv(file,
                                  col_types = readr::cols(
                                      .default = col_integer(),
                                      X1 = col_character()
                                  ))
        counts <- dplyr::select(counts, -X1)
        counts <- as.matrix(counts)
    } else {
        stop("Dataset not valid")
    }

    return(counts)
}
