\name{r2.nakagawa.glm}
\alias{r2.nakagawa.glm}

\title{
Goodness-of-fit measure for binary/logistic model
}
\description{
GOF measure for GLMs are generally pretty bad.  For binary data, likelihood based methods will be low valued.  Nakagawa or MZ type goodness-of-fit measures are based on the latent continuous probability that you're trying to predict.
}
\usage{
r2.nakagawa.glm(m, varType = 2)
}
\arguments{
  \item{m}{
  GLM or GLMM
}
  \item{varType}{
  Type of variance to use, could be a single number or a vector. (see Details)
}
}
\details{
This GOF measure is based on variance based on different parts of the model.  Pretty straight forward if it's a GLM, but if it's a GLMM, there are several different options of variances to use that'll be slightly different.  Both marginal and conditional versions are returned, but VarR could either come from the generalized variance (varType=1), the total variance (varType=2), or based on the entire correlation structure of fixed with randome effects (i.e. var(fitted(m)); varType=3). 

Make sure you check if there's a weird random effect structure.  That'll need to be dealt with separately.  
}
\value{
R2 value for the model 
If class(m)=='mer', output will be a list with marginal and conditional R2
}
\references{
Nakagawa and Schielzeth. 2013. "A general and simple method for obtaining R2 from generalized linear mixed-effects models." Methods in Ecology and Evolution.
}
\author{
ME
}

