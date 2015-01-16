\name{file.clean}
\alias{file.clean}

\title{
Clean out extra files from a knitr/latex build
}
\description{
Clean out extra files from a knitr/latex build
}
\usage{
file.clean(fname)
}

\arguments{
  \item{fname}{
base file name for the knitr/tex files
}
}
\details{
Removes everything with fname except the .rnw, .tex, and .pdf.
}
\value{
Returns file names that got dropped.
}

\author{
ME
}
