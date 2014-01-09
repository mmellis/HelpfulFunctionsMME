\name{tryN}
\alias{tryN}

\title{
Adapted try() function
}
\description{
  A shortened version of tryCatch that will return NA on error
}
\usage{
tryN(f)
}

\arguments{
  \item{f}{
Something you want to evaluate that might fail
}
}
\details{
tryCatch(f, error = function(e) return(NA))
}
\value{
NA if f returns an error
}

