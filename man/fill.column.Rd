\name{fill.column}
\alias{fill.column}

\title{
Fill down column values
}
\description{
 In Excel, when you don't want to repeat column values over and over again in block.  This will copy values down until it hits the next known value.
}
\usage{
fill.column(x)
}
\arguments{
  \item{x}{
 a vector with a lot of NAs
}
}
\details{
  Whatever is the last non-NA value will be copied until it hits another non-NA value.
}
\value{
 a new vector with no NAs
}

\author{
ME
}

\examples{
  x<-c(1, NA, NA, NA, 2, NA ,NA, 3, NA,NA,NA, 4,5,3, NA,NA)
  cbind(x, fill.column(x))
  }
