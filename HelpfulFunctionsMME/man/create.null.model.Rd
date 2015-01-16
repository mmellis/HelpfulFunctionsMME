\name{create.null.model}
\alias{create.null.model}

\title{
Create a null model
}
\description{
Create a null model with just an intercept term
}
\usage{
create.null.model(fit)
}

\arguments{
  \item{fit}{
A model
}
}

\value{
The null model for that model
}

\author{
ME
}

\examples{
#function (fit) 
#{
#    form <- paste(formula(fit))
#    if (class(fit)[1] == "mer") 
#        form[3] <- gsub(".*\\+ \\(", "1 + (", form[3])
#    else form[3] <- 1
#    form <- formula(paste(form[2], form[3], sep = "~"))
#    fit <- update(fit, form)
#    return(fit)
#  }
}
