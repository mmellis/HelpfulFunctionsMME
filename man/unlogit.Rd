\name{unlogit}
\alias{unlogit}

\title{
Backtransforms logit data
}
\description{
Backtransforms logit data
}
\usage{
unlogit(x)
}
\details{
exp(x)/(1+exp(x))
}

\arguments{
  \item{x}{
a numeric vector
}
}

\value{
backtransformed values
}

\author{
ME
}

\examples{
x<- -100:100
unlogit(x)
}