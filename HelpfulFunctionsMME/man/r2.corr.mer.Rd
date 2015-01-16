\name{r2.corr.mer}
\alias{r2.corr.mer}

\title{
Coefficient of determiniation for (G)LMMs
}
\description{
Goodness of fit for (G)LMM based on individual effects - Conditional
}
\usage{
r2.corr.mer(m)
}

\arguments{
  \item{m}{
 (G)LMM object 
}
}
\details{
Correlation between responses and fitted values.  Conditional = based on individual random effects.

From Doug Bates: "Assuming that one wants to define an R^2 measure, I think an
argument could be made for treating the penalized residual sum of
squares from a linear mixed model in the same way that we consider the
residual sum of squares from a linear model. Or one could use just
the residual sum of squares without the penalty or the minimum
residual sum of squares obtainable from a given set of terms, which
corresponds to an infinite precision matrix. I don't know, really.
It depends on what you are trying to characterize."
}
\value{
R2
}
\references{
glmm.wikidot.com/faq
}
\author{
Written by Jarrett Byrnes
}
