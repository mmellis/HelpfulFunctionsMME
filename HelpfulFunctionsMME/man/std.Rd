\name{std}
\alias{std}

\title{
Standardize a variable by its mean and sd
}
\description{ Standardize a variable by its mean and sd}

\usage{
std(x)
}
\arguments{
  \item{x}{
 a numeric vector
}
}

\value{
Same thing, standardized by its mean and sd.  NAs removed.
}

\author{
ME
}
\examples{
x<-rnorm(100,mean=10, sd=12)
std(x)
}