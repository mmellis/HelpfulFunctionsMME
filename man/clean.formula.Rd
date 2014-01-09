\name{clean.formula}
\alias{clean.formula}
\title{
  Makes a shorter version of a formula for display
}
\description{
Used in stepAICc and others.  If you want to use formulas from a big AIC analysis in a table, they need to be shorter and prettier.  This sorts the terms, replaces the top model form with '...', etc.  
}
\usage{
clean.formula(x, left = F, base.form = NULL)
}

\arguments{
  \item{x}{
    A formula, either as a character string or a factor.  Or it could be a vector of formulas.
}
  \item{left}{
   A logical, do you want the most common/consistent terms on the left side of the formula or the right.  How do you want the formulas to line up?
}
  \item{base.form}{
  I don't think this does anything anymore
}
}
\details{
  Takes the formula.  Turns it into text.  Substitutes out spaces and removes the response variable.  Splits the strings by terms ('+').  Finds interaction terms and drops repeated terms, switch to '*' notation instead of including both interaction and main effects in formula.  Sorts the terms by their frequency of appearance in the vector of formulas. Sort and return the cleaned up vector of formulas.
}
\value{
  A much neater formula
}
\author{
Martha
}
\examples{
x<-c("SURV ~ Denali.snow + Denali.snow:initial.weight + initial.weight + summer.rain.yr1", 
"SURV ~ Denali.snow + Denali.snow:initial.weight + initial.weight", 
"SURV ~ Denali.snow + Denali.snow:initial.weight + initial.weight + summer.rain.yr2"
)
clean.formula(x)
}
