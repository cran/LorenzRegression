% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LR-coef.R
\name{coef.LR}
\alias{coef.LR}
\title{Estimated coefficients for the Lorenz Regression}
\usage{
\method{coef}{LR}(object, ...)
}
\arguments{
\item{object}{Output of a call to \code{\link{Lorenz.Reg}}, where \code{penalty="none"}.}

\item{...}{Additional arguments.}
}
\value{
a vector gathering the estimated coefficients
}
\description{
\code{coef.LR} provides the estimated coefficients for an object of class \code{LR}.
}
\examples{
data(Data.Incomes)
NPLR <- Lorenz.Reg(Income ~ ., data = Data.Incomes, penalty = "none")
coef(NPLR)

}
\seealso{
\code{\link{Lorenz.Reg}}
}
