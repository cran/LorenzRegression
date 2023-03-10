% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLR-print.R
\name{print.PLR}
\alias{print.PLR}
\title{Printing method for the Penalized Lorenz Regression}
\usage{
\method{print}{PLR}(x, ...)
}
\arguments{
\item{x}{Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty!="none"}.}

\item{...}{Additional arguments.}
}
\value{
No return value, called for printing an object of class \code{PLR} to the console
}
\description{
\code{print.PLR} prints the arguments and estimated coefficients of an object of class \code{PLR}.
}
\examples{
data(Data.Incomes)
PLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "SCAD",
                  sel.choice = c("BIC","CV"), h.grid = nrow(Data.Incomes)^(-1/5.5),
                  eps = 0.01, seed.CV = 123, nfolds = 5)
print(PLR)

}
\seealso{
\code{\link{Lorenz.Reg}}
}
