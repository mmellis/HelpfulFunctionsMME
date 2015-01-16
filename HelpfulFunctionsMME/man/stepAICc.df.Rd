\name{stepAICc.df}
\alias{stepAICc.df}

\title{
stepAICc.df
}
\description{
Produce a dataframe for a set of models
}
\usage{
stepAICc.df(fits)
}

\arguments{
  \item{fits}{
list of output from stepAICc()
}
}

\value{
data.frame
}

\author{
ME
}

\examples{
#function (fits) 
#{
#    fit.df <- rbind.fill(lapply(fits, function(x) as.data.frame(x[1:5])))
#    null.deviance <- deviance(update(fits[[1]]$object, ".~1"))
#    fit.df$r2.verhoef <- sapply(fits, function(x) {
#        1 - x$deviance/null.deviance
#    })
#    fit.df$r2.nakagawa <- sapply(fits, function(x) {
#        if (!is.null(x$object)) {
#            r2.nakagawa.glm(x$object)
#        }
#        else {
#            r2.nakagawa.glm(update(fits[[1]]$object, x$formula))
#        }
#    })
#    fit.df$k <- sapply(strsplit(paste(fit.df$formula), split = "+", 
#        fixed = T), function(x) {
#        x <- gsub(".*~", "", x)
#        x <- gsub(" ", "", x)
#        if (all(x != "1")) 
#            x <- c(x, "1")
#        length(x)
#    })
#    fit.df$formula <- clean.formula(fit.df$formula)
#    fit.df$dAIC <- fit.df$AIC - fit.df$AIC[1]
#    fit.df$ModelID <- 1:nrow(fit.df)
#    fit.df <- subset(fit.df, select = c(deviance, df.resid, k, 
#        AIC, dAIC, r2.verhoef, r2.nakagawa, formula, ModelID))
#    return(fit.df)
#  }
}