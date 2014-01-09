\name{substr1}
\alias{substr1}

\title{
Substr() of just one element
}
\description{
Substr() of just one element
}
\usage{
substr1(x, ch = 1)
}

\arguments{
  \item{x}{
  a character vector
}
  \item{ch}{
  which element of the string.  ch<0 that element will get dropped.
}
}

\value{
If ch>0, the ch-th character of x.
If ch<0, a substr of x without the ch-th character.
}

\author{
ME
}
\examples{
substr1('hello',3)
substr1('hello',-3)
substr1('hello',-c(1:3))
substr1(c('hello','goodbye'),3)
}