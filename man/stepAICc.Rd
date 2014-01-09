\name{stepAICc}
\alias{stepAICc}

\title{
Stepwise AIC
}
\description{
Build a set of candidate models for a stepwise AIC model selection
}
\usage{
stepAICc(object, scope, scale = 0, direction = c("both", "backward", "forward"), trace = 1, steps = 1000, use.start = FALSE, k = 2, bandwidth = 10, marginal = NULL, AICc = TRUE, delta = 2, delta2 = Inf, joint = NULL, positive.definite = FALSE, ...)
}

\arguments{
  \item{object}{
Initial model
}
  \item{scope}{
list(lower=formula(low), high=formula(global model))
}
  \item{scale}{
??
}
  \item{direction}{
"both", "backward", or "forward"
}
  \item{trace}{
How much info to output.
}
  \item{steps}{
Limit on how much to search
}
  \item{use.start}{
??
}
  \item{k}{
parameters to use in AIC, I think?
}
  \item{bandwidth}{
??
}
  \item{marginal}{
Has to do with you to search.  Leave blank.  Don't break rules for including interactions without main effects.
}
  \item{AICc}{
T/F use AICc instead of AIC
}
  \item{delta}{
How big of a jump is needed to flag a model to expand
}
  \item{delta2}{
Which models to output
}
  \item{joint}{
??
}
  \item{positive.definite}{
??
}
  \item{\dots}{
other stuff??
}
}

\value{
list of output with formulas for models it tested and often times other information
}

\author{
ME (edited a lot from stuff found on the internet)
}
