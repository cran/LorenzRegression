#' Plots for the Lorenz regression
#'
#' \code{autoplot} generates a plot for an object of class \code{"LR"} and returns it as a \code{ggplot} object.
#' The \code{plot} method is a wrapper around \code{autoplot} that directly displays the plot,
#' providing a more familiar interface for users accustomed to base R plotting.
#'
#' @aliases plot.LR autoplot.LR_boot plot.LR_boot
#' @param x An object of class \code{"LR"}.
#' @param object An object of class \code{"LR"}.
#' @param ... Additional arguments passed to \code{\link{Lorenz.graphs}}.
#'
#' @return \code{autoplot} returns a \code{ggplot} object representing the Lorenz curve of the response and the concentration curve of the response with respect to the estimated index. \code{plot} directly displays this plot.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' ## For examples see example(Lorenz.Reg)
#'
#' @importFrom ggplot2 ggtitle autoplot
#' @importFrom stats update.formula
#'
#' @method autoplot LR
#' @export

autoplot.LR <- function(object, ...){

  if (is.null(object$theta)) stop("No plots are available for an empty model.")

  formula <- update.formula(object, . ~ index)
  data <- data.frame(object$y,predict.LR(object))
  names(data) <- all.vars(formula)

  g <- Lorenz.graphs(formula, data, weights = object$weights, ...)
  g <- g + ggtitle("Observed and explained inequality")

  g

}

#' @importFrom graphics plot
#' @method plot LR
#' @rdname autoplot.LR
#' @export
plot.LR <- function(x, ...) {
  print(autoplot(x, ...))
}

#' @method autoplot LR_boot
#' @export
autoplot.LR_boot <- function(object, ...){
  NextMethod("autoplot")
}

#' @method plot LR_boot
#' @export
plot.LR_boot <- function(x, ...) {
  print(autoplot(x, ...))
}

