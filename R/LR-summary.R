#' Summary for the Lorenz Regression
#'
#' \code{summary.LR} provides a summary for an object of class \code{LR}.
#'
#' @param object Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty=="none"}.
#' @param ... Additional arguments
#'
#' @return A summary displaying the explained Gini coefficient, Lorenz-\eqn{R^2} and a table gathering the estimated coefficients, including p-values if bootstrap was performed.
#'
#' @seealso \code{\link{Lorenz.Reg}}
#'
#' @examples
#' data(Data.Incomes)
#' NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
#' summary(NPLR)
#'
#' @import knitr
#'
#' @method summary LR
#' @export

summary.LR <- function(object, ...){

  LR <- object
  theta.mat <- as.matrix(LR$theta)
  if ("pval.theta" %in% names(LR)){
    theta.mat <- cbind(theta.mat, LR$pval.theta)
    txt <- "Estimated coefficients and associated p-values"
    colnames(theta.mat) <- c("estimate","p-value")
  }else{
    txt <- "Estimated coefficients"
    colnames(theta.mat) <- c("estimate")
  }

  theta.table <- knitr::kable(theta.mat)

  cat(paste0("The explained Gini coefficient is of ",round(LR$Gi.expl,5)),
            "",
            paste0("The Lorenz-R2 is of ",round(LR$LR2,5)),
            "",
            txt,
            theta.table,
            sep = "\n")

}
