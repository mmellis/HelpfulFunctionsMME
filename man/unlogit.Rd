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

Using a straight backtransformation for the logit function can lead to bias if applied to a distribution summary, such as a mean or CI limit. This bias is assumed to be small, so technically we're really apply a 1st order Delta method correction.  Tedious, but sometimes the bias can be important and it may be necessary to use a higher order taylor series expansion.
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
