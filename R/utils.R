#' Character round
#'
#' Round number for presentation. Keeps the same number of digits regardless
#' of what they are, so 0.101 will be rounded to 0.10 not 0.1 like the regular
#' round function.
#'
#' @param x Vector of numbers to round
#' @param digits Number of digits to round to
#'
#' @return Rounded character vector
chrRound <- function(x, digits = 1) {

    if(digits < 1) {
        stop("This is intended for the case digits >= 1.")
    }

    if(length(digits) > 1) {
        digits <- digits[1]
        warning("Digits should be a single integer. Only using digits[1].")
    }

    rounded <- sprintf(paste("%.", digits, "f", sep = ""), x)

    # Convert "-0.00" to "0.00"
    zero <- paste0("0.", paste(rep("0", digits), collapse = ""))
    rounded[rounded == paste0("-", zero)] <- zero

    return(rounded)
}


#' Logistic function
#'
#' Implementation of the logistic function
#'
#' @param x value to apply the function to.
#' @param x0 midpoint parameter. Gives the centre of the function.
#' @param k shape parameter. Gives the slope of the function.
#'
#' @return Value of logistic funciton with given parameters
logistic <- function(x, x0, k) {
    1 / (1 + exp(-k * (x - x0)))
}

#' Complete list of MCRI palettes
#'
#' Use \code{\link{mcriPalette}} to construct palettes of desired length.
#'
#' @export
mcri.palettes <- list(
    themes = c("#EC008C", "#00ADEF", "#8DC63F", "#00B7C6", "#F47920",
               "#7A52C7"),
    themesMid = c("#F499C2", "#6BCFF6", "#C4DF9B", "#92D6DE", "#FAB783",
                  "#B09ECB"),
    themesLite = c("#F9CCDF", "#B9E5FA", "#DFEDCB", "#C8E8ED", "#FDD8BB",
                   "#D3CAE3"),
    themesPaired = c("#EC008C", "#F499C2", "#00ADEF", "#6BCFF6", "#8DC63F",
                     "#C4DF9B", "#00B7C6", "#92D6DE", "#F47920", "#FAB783",
                     "#7A52C7", "#B09ECB"),
    themesTripled = c("#EC008C", "#F499C2", "#F9CCDF", "#00ADEF", "#6BCFF6",
                      "#B9E5FA", "#8DC63F", "#C4DF9B", "#DFEDCB", "#00B7C6",
                      "#92D6DE", "#C8E8ED", "#F47920", "#FAB783", "#FDD8BB",
                      "#7A52C7", "#B09ECB", "#D3CAE3"),
    blues = c("#00A5D2", "#0053A1", "#092F5E"),
    bluesMid = c("#6FCEE7", "#78A7D7", "#717BAB"),
    bluesLite = c("#BAE4F2", "#B7CDE8", "#ACB1D0"),
    bluesPaired = c("#00A5D2", "#6FCEE7", "#0053A1", "#78A7D7", "#092F5E",
                    "#717BAB"),
    bluesTripled = c("#00A5D2", "#6FCEE7", "#BAE4F2", "#0053A1", "#78A7D7",
                     "#B7CDE8", "#092F5E", "#717BAB", "#ACB1D0"),
    logo = c("#092F5E", "#0053A1", "#00A5D2", "#EC008C")
)

#' A MCRI palette generator
#'
#' These are a handful of colour palettes based on the MCRI branding.
#'
#' @param n Number of colors desired. If omitted, uses all colours.
#' @param name Name of desired palette. Choices are:
#'   \code{MCRI}
#' @param type Either "continuous" or "discrete". Use continuous if you want
#'   to automatically interpolate between colours.
#' @return A vector of colours.
#'
#' @export
#' @keywords colours, colors
#' @examples
#' mcriPalette("MCRI")
#'
#' # If you need more colours than normally found in a palette, you
#' # can use a continuous palette to interpolate between existing
#' # colours
#' pal <- mcriPalette(21, name = "MCRI", type = "continuous")
#' image(volcano, col = pal)
mcriPalette <- function(name, n, type = c("discrete", "continuous")) {
    type <- match.arg(type)

    pal <- mcri.palettes[[name]]
    if (is.null(pal))
        stop("Palette not found.")

    if (missing(n)) {
        n <- length(pal)
    }

    if (type == "discrete" && n > length(pal)) {
        stop("Number of requested colours greater than what palette can offer")
    }

    out <- switch(type,
                  continuous = colorRampPalette(pal)(n),
                  discrete = pal[1:n]
    )
    structure(out, class = "palette", name = name)
}
