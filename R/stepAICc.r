stepAICc <- function (object, scope, scale = 0, direction = c("both", 
"backward", "forward"), trace = 1, steps = 1000, use.start = FALSE, k = 
2, bandwidth = 10, marginal = NULL, AICc = TRUE, delta = 2 ,delta2 = Inf, joint = 
NULL, positive.definite = FALSE, ...){ 
### Subfunctions    
    mydeviance <- function(x, ...) { 
        dev <- deviance(x) 
        if (!is.null(dev)){ 
            dev 
        } else { 
            extractAIC(x, k = 0)[2] 
        } 
    } 
#    
    cut.string <- function(string) { 
        if (length(string) > 1){ 
            string[-1] <- paste("\n", string[-1], sep = "") 
        } 
        string 
    } 
#    
    rename.interactions<-function(ff, sepM=":"){     #Renames interactions in alpha-order (B:A => A:B)
      ff<-strsplit(ff, sepM,fixed=T)
      sapply(ff, function(x) paste(sort(x), collapse = sepM))
    }
#    
    sort.terms<-function(terms.factor){
      if(class(terms.factor)=='matrix'){
              colnames(terms.factor) <- rename.interactions(colnames(terms.factor)) 
              terms.factor <- terms.factor[, order(colnames(terms.factor))] 
            } else { 
              terms.factor <- numeric(0) 
            } 
      return(terms.factor)
    } 
#    
    sort.fixed<-function(fixed){
        fixed<-gsub(' ','',fixed)
       fixed<-strsplit(fixed, '+', fixed=T)
       fixed<-lapply(fixed, function(x) sort(rename.interactions(x)))
       unlist(lapply(fixed, paste, collapse=" + "))
       }
#           
    text.formula<-function(fit){
      use.form<-formula(fit)
      use.form<-paste(use.form)
      use.form[3]<-sort.fixed(use.form[3])
      paste(use.form[c(2,1,3)],collapse=' ')
    } 
#    
    factor.scope <- function (factor, scope, marginal = NULL, joint = NULL) 
{ 
    drop <- scope$drop 
    add <- scope$add 
    if (length(factor) && !is.null(drop)) { 
        facs <- factor 
        if (length(drop)) { 
            if(is.null(nrow(drop))){nmdrop <- names(drop)[-1]}else{nmdrop <- colnames(drop)} 
            nmfac <- colnames(factor) 
            nmfac0 <- sort.fixed(nmfac)
            nmdrop0 <- sort.fixed(nmdrop)
            where <- match(nmdrop0, nmfac0, 0) 
            if (any(!where)){ 
                stop(gettextf("lower scope has term(s) %s not included in model",  paste(sQuote(nmdrop[where == 0]), collapse = ", ")), domain = NA) 
            } 
            facs <- factor[, -where, drop = FALSE] 
            nmdrop <- nmfac[-where] 
        } else { 
            nmdrop <- colnames(factor) 
        } 
        if (ncol(facs) > 1) { 
            keep <- rep.int(TRUE, ncol(facs)) 
            f <- crossprod(facs > 0) 
            for (i in seq(keep)){ 
                keep[i] <- max(f[i, -i]) != f[i, i] 
            } 
            nmdrop <- nmdrop[keep] 
        } 
        keep <- unique(unlist(lapply(strsplit(nmdrop, ":", fixed = TRUE), function(x){ 
            whereMarginal <- match(x, names(marginal)) 
            if(!all(is.na(whereMarginal))){ 
                if(length(x) == 1){ 
                    marginal[[whereMarginal]] 
                } else { 
                    lowerOrderTerms <- lapply(seq_along(whereMarginal), 
function(z){ 
                        if(is.na(whereMarginal[z])){ 
                           x[z] 
                        } else { 
                            c(x[z], marginal[[whereMarginal[z]]]) 
                        } 
                    }) 
                    lowerOrder <- apply(expand.grid(lowerOrderTerms), 1, 
function(z){ 
                        paste(sort(z), collapse = ":") 
                    }) 
                    lowerOrder[!lowerOrder %in% paste(sort(x), collapse = ":")] 
                } 
            } 
        }))) 
        nmdrop <- nmdrop[!nmdrop %in% keep] 
        if(length(joint) & length(nmdrop)){ 
            nmdrop <- sapply(strsplit(nmdrop, ":", fixed = TRUE), 
function(sdrop){ 
                whereJoint <- sapply(joint, function(x){ 
                    !all(!sdrop %in% x) 
                }) 
                if(sum(whereJoint) > 0){ 
                    if(length(sdrop) == 1){ 
                        paste(unlist(joint[whereJoint]), collapse = "-") 
                    } else { 
                        whichJoint <- !sdrop %in% unlist(joint[whereJoint]) 
                        paste(sdrop[whichJoint], unlist(joint[whereJoint]), sep = ":", collapse = "-") 
                    } 
                } else { 
                    paste(sdrop, collapse = ":") 
                } 
            }) 
            nmdrop <- unique(nmdrop) 
        } 
    } else { 
        nmdrop <- character(0) 
    } 
    if (!length(add)){ 
        nmadd <- character(0) 
    } else { 
        nmfac <- colnames(factor) 
        nmadd <- colnames(add) 
        if (!is.null(nmfac)) { 
            nmfac0 <- sort.fixed(nmfac) 
            nmadd0 <- sort.fixed(nmadd) 
            where <- match(nmfac0, nmadd0, 0) 
            if (any(!where)){ 
                stop(gettextf("upper scope does not include model term(s) %s", paste(sQuote(nmfac[where == 0]), collapse = ", ")), domain = NA) 
            } 
            nmadd <- nmadd[-where] 
            add <- add[, -where, drop = FALSE] 
        } 
        if (ncol(add) > 1) { 
            keep <- rep.int(TRUE, ncol(add)) 
            f <- crossprod(add > 0) 
            for (i in seq(keep)){ 
                keep[-i] <- keep[-i] & (f[i, -i] < f[i, i]) 
            } 
            nmadd <- nmadd[keep] 
        } 
        if(length(nmadd)){ 
            keep <- sapply(strsplit(nmadd, ":", fixed = TRUE), 
function(x){ 
                whereMarginal <- match(x, names(marginal)) 
                if(!all(is.na(whereMarginal))){ 
                    if(length(x) == 1){ 
                        all(!marginal[[whereMarginal]] %in% nmadd) 
                    } else { 
                        lowerOrderTerms <- lapply(seq_along(whereMarginal), function(z){ 
                            if(is.na(whereMarginal[z])){ 
                               x[z] 
                            } else { 
                                c(x[z], marginal[[whereMarginal[z]]]) 
                            } 
                        }) 
                        lowerOrder <- apply(expand.grid(lowerOrderTerms), 1, function(z){ 
                            paste(sort(z), collapse = ":") 
                        }) 
                        lowerOrder <- lowerOrder[!lowerOrder %in% paste(sort(x), collapse = ":")] 
                        all(!lowerOrder %in% nmadd) 
                    } 
                } else { 
                    TRUE 
                } 
            }) 
            nmadd <- nmadd[keep] 
        } 
        if(length(nmadd) & length(joint)){ 
            nmadd <- sapply(strsplit(nmadd, ":", fixed = TRUE), 
function(sadd){ 
                whereJoint <- sapply(joint, function(x){ 
                    !all(!sadd %in% x) 
                }) 
                if(sum(whereJoint) > 0){ 
                    whichJoint <- !sadd %in% unlist(joint[whereJoint]) 
                    if(length(sadd) == 1){ 
                        paste(unlist(joint[whereJoint]), collapse = "+") 
                    } else { 
                        paste(sadd[whichJoint], 
unlist(joint[whereJoint]), sep = ":", collapse = "+") 
                    } 
                } else { 
                    paste(sadd, collapse = ":") 
                } 
            }) 
            nmadd <- unique(nmadd) 
            
        } 
    } 
    list(drop = nmdrop, add = nmadd) 
} 
### End Subfunctions
  
    Terms <- terms(object) 
    object$formula <- Terms 
    if (inherits(object, "lme")){ 
        object$call$fixed <- Terms 
    } else { 
        if (inherits(object, "gls")){ 
            object$call$model <- Terms 
        } else { 
            object$call$formula <- Terms 
        } 
    } 
    if (use.start){ 
        warning("'use.start' cannot be used with R's version of glm") 
    } 
    md <- missing(direction) 
    direction <- match.arg(direction) 
    backward <- direction == "both" | direction == "backward" 
    forward <- direction == "both" | direction == "forward" 
    if (missing(scope)) { 
        fdrop <- numeric(0) 
        fadd <- attr(Terms, "factors") 
        if (md){ forward <- FALSE } 
    } else { 
        if (is.list(scope)) { 
            if(!is.null(fdrop <- scope$lower)){ 
                fdrop <- attr(terms(update.formula(object, fdrop)), "factors") 
                fdrop <- sort.terms(fdrop) 
            }
            if(!is.null(fadd <- scope$upper)){ 
                fadd <- attr(terms(update.formula(object, fadd)), "factors") 
                fadd <- sort.terms(fadd)
            } 
        } else { 
            if(!is.null(fadd <- scope)){ 
                fadd <- attr(terms(update.formula(object, scope)), "factors") 
                fadd <- sort.terms(fadd) 
            } 
            fdrop <- numeric(0) 
        } 
    } 
    
    if(length(fdrop)==0) object<-update(object, formula=~1)
    
    models <- vector("list", steps) 
    if (is.list(object) && (nmm <- match("nobs", names(object), 0)) > 0) { 
        n <- object[[nmm]] 
    } else { 
        n <- length(residuals(object)) 
    } 
    fit <- object 
        bAIC <- extractAIC(fit, scale, k = k, ...) 
#	bAIC <- extractAIC(fit, scale, k = k) 
        edf <- bAIC[1] 
    if(AICc){ 
        bAIC <- bAIC[2] + (2 * edf + (edf * 1)) / (n - edf - 1) 
    } else { 
        bAIC <- bAIC[2] 
    } 
    if (is.na(bAIC)){ 
        stop("AIC is not defined for this model, so stepAIC cannot proceed") 
    } 
    Terms <- terms(fit) 
# reorder the covariates in alphabetical order 
    response <- paste(formula(fit))[2]
    fixed <- sort.fixed(paste(formula(fit))[3])
     
    nm <- 1 
#append the basic model to the list of the results
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - edf, AIC = bAIC, formula = text.formula(fit), expanded = FALSE, object = fit) 
    bestAIC <- bAIC 
    selection <- 1 
    while (nm < steps & length(selection) > 0) { 
#create a list of formulas which have to be expanded 
        formulas <- sapply(models[selection], function(x){x$formula}) 
#create a vector with AIC values of the model to be expanded. 
        cAIC <- sapply(models[selection], function(x){x$AIC}) 
        if(trace >= 1){ 
            cat("\n") 
            cat(nm, "models tested.\n") 
            if(length(cAIC) > 1){ 
                cat(length(cAIC), "models to expand. AIC range:", round(range(cAIC[is.finite(cAIC)], na.rm = TRUE), 3), "\n\n") 
            } else { 
                cat("1 model to expand.\n\n") 
            } 
        } 
#set the expanded flag on the models selected to expand 
        models[selection] <- lapply(models[selection], function(x){ 
            x$expanded <- TRUE 
            if(is.infinite(x$AIC)){ 
                x$AIC <- NA 
            } 
            x 
        }) 
#expand each selected model 
        for(i in seq_along(formulas)){ 
            form <- formulas[i] 
            if(trace >= 2){ 
                cat("Expand model", i, "\n", 
cut.string(deparse(as.vector(formula(form)))), "\nAIC:", round(cAIC[[i]], 3), " Best AIC:", round(bestAIC, 3), "AIC range:", round(range(cAIC[is.finite(cAIC)], na.rm = TRUE), 3), "\n\n") 
                utils::flush.console() 
            } 
#select the terms which can be added of dropped 
            Terms <- terms(as.formula(form)) 
            ffac <- attr(Terms, "factors") 
            if (!is.null(sp <- attr(Terms, "specials")) && !is.null(st <- sp$strata)){ 
                ffac <- ffac[-st, ] 
            } 
            scope <- factor.scope(ffac, list(add = fadd, drop = fdrop), marginal = marginal, joint = joint) 
#to backward selection
            if (backward && length(scope$drop)) { 
                  rhs<-sapply(scope$drop, function(x) update.formula(form, formula(paste('. ~ . ', x, sep='-'))))
                    rhs<-sapply(rhs, function(rhst) formula(text.formula(rhst)))
                  rhs.drop<-sapply(rhs, function(rhst){ all(sapply(models[seq_len(nm)], function(x){formula(x$formula) != rhst})) })
#ignor smaller models which have allready been calculated 
                scope$drop <- scope$drop[rhs.drop] 
#drop all suitable terms one by one 
                if(all(!is.na(scope$drop)) & length(scope$drop)){ 
                    for(sdrop in scope$drop){ 
                        newForm <- update.formula(form, as.formula(paste(". ~ . ", sdrop, sep = "-")))
                        newObject <- update(fit, newForm, evaluate = FALSE) 
                        newObject <- try(eval.parent(newObject)) 
#check if the model gave an error or not 
                        if(class(newObject)[1] != "try-error"){ 
                            bAIC <- extractAIC(newObject, scale, k = k, ...) 
#                            bAIC <- extractAIC(newObject, scale, k = k) 
                            edf <- bAIC[1] 
                            if(AICc){ 
                                bAIC <- bAIC[2] + (2 * edf + (edf * 1)) / (n - edf - 1) 
                            } else { 
                                bAIC <- bAIC[2] 
                            } 
                            deviance <- mydeviance(newObject) 
                            if(positive.definite && class(try(intervals(newObject))) == "try-error"){ 
                                bAIC <- ifelse(is.infinite(cAIC[i]), NA, -Inf) 
                            } else {                           
                                if(bestAIC > bAIC){ 
                                    bestAIC <- bAIC 
                                } 
                            } 
                        } else { 
#set the AIC temporary to -Inf if the model gave an error. This will be 
#reset later on to NA. It is set to -Inf so the model will have the 
#lowest AIC and will be expanded. 
                            edf <- n 
                            bAIC <- ifelse(is.infinite(cAIC[i]), NA, -Inf) 
                            newObject <- NULL 
                            deviance <- NA 
                        } 
                        newFixed<-sort.fixed(paste(newForm)[3])
                        if(class(attr(terms(newForm), "factors"))[1]=='matrix'){  
                            expanded <- FALSE 
                        } else {  
                            expanded <- TRUE 
                            bAIC <- ifelse(is.infinite(bAIC), NA, bAIC) 
                        } 
                        nm <- nm + 1 
                        models[[nm]] <- list(deviance = deviance, df.resid = n - edf, AIC = bAIC, formula = text.formula(formula(newObject)), expanded = expanded, object = newObject) 
                        if(trace >= 3){ 
                            cat("Model", nm, "Drop: ", sdrop, "\n", 
cut.string(deparse(as.vector(as.formula(paste(response, newFixed, sep = " ~ "))))), "\nAIC: ", round(bAIC, 3), " Best AIC: ", round(bestAIC, 3), "\n\n") 
                            utils::flush.console() 
                        } 
                    } 
                } 
            }
            if (forward && length(scope$add)) {                 
                rhs<-sapply(scope$add, function(x) update.formula(form, formula(paste('. ~ . ', x, sep='+'))))
                    rhs<-sapply(rhs, function(rhst) formula(text.formula(rhst)))
                rhs.drop<-sapply(rhs, function(rhst){ all(sapply(models[seq_len(nm)], function(x){formula(x$formula) != rhst}))})
                scope$add <- scope$add[rhs.drop] 
                if(all(!is.na(scope$add)) & length(scope$add)){ 
                    for(sadd in scope$add){ 
                        newForm <- update.formula(form, as.formula(paste(". ~ . ", sadd, sep = "+"))) 
                        newObject <- update(fit, newForm, evaluate = FALSE) 
                        newObject <- try(eval.parent(newObject)) 
                        if(class(newObject)[1] != "try-error"){ 
                            bAIC <- extractAIC(newObject, scale, k = k, ...) 
#                            bAIC <- extractAIC(newObject, scale, k = k) 
                            edf <- bAIC[1] 
                            if(AICc){ 
                                bAIC <- bAIC[2] + (2 * edf + (edf * 1)) / (n - edf - 1) 
                            } else { 
                                bAIC <- bAIC[2] 
                            } 
                            deviance <- mydeviance(newObject) 
                            if(positive.definite && class(try(intervals(newObject))) == "try-error"){ 
                                bAIC <- ifelse(is.infinite(cAIC[i]), NA, -Inf) 
                            } else {                           
                                if(bestAIC > bAIC){ 
                                    bestAIC <- bAIC 
                                } 
                            } 
                        } else { 
                            edf <- NA 
                            deviance <- NA 
                            bAIC <- ifelse(is.infinite(cAIC[i]), NA, -Inf) 
                            newObject <- NULL 
                        } 
                        newFixed <- sort.fixed(paste(newForm)[3])
                        nm <- nm + 1 
                        models[[nm]] <- list(deviance = deviance, df.resid = n - edf, AIC = bAIC, 
                            formula = text.formula(formula(newObject)), expanded = FALSE, object = newObject) 
                        if(trace >= 3){ 
                            cat("Model", nm, "Add: ", sadd, "\n", 
                              cut.string(deparse(as.vector(as.formula(paste(response, newFixed, sep = " ~ "))))), 
                                "\nAIC: ", round(bAIC, 3), " Best AIC: ", round(bestAIC, 3), "\n\n") 
                            utils::flush.console() 
                        } 
                    } 
                } 
            } 
        } 
        models[seq_len(nm)] <- models[seq_len(nm)][order(sapply(models[seq_len(nm)], 
function(x){x$AIC}))] 
        selection <- which(sapply(models[seq_len(nm)], function(x){ 
            (x$AIC - bestAIC < delta) | is.infinite(x$AIC) 
        }) == TRUE) 
        if(length(selection) < bandwidth){ 
            selection <- seq_len(min(c(nm, bandwidth))) 
        } 
        unselect <- seq_len(nm)[!seq_len(nm) %in% selection] 
        models[unselect] <- lapply(models[unselect], function(x){ 
            x$object <- NULL 
            x 
        }) 
        selection <- selection[sapply(models[selection], 
function(x){!x$expanded})] 
    } 
        if(nm > steps){ 
                cat("Warning: algorithm did not converge") 
                models[[1]]$Converged <- FALSE 
        } else { 
                models[[1]]$Converged <- TRUE 
        } 

    models<-models[sapply(models, function(x){!is.null(x)})]    
    models<-models[(sapply(models, function(x) x$AIC)-models[[1]]$AIC)<=delta2] 
    return(models)
} 
##############################     End Main Function

#######################

  clean.formula<-function(x,left=F, base.form=NULL){
    # Collapses interactions
    if(class(x)=='factor') x<-paste(x)
    
        x<-gsub(' ','',x)
        x<-gsub('.*~','',x)
        
        ls.x<-strsplit(x,'\\+')
        ls.x<-lapply(ls.x, function(xx) {
            spl<-grep(':',xx)
            if(length(spl)>0){
              int.terms<-unique(unlist(strsplit(xx[spl], ':')))
              xx<-xx[-match(int.terms, xx)]  #Removes main effect terms from disply.
              xx<-gsub(':','*',xx)
              }
            return(xx)})
        
        # Sort terms
        all.terms<-as.data.frame(table(unlist(ls.x)),stringsAsFactors=F)
          if(left){ 
            use.order<-order(-all.terms$Freq, nchar(all.terms$Var1))
          } else {use.order<-order(all.terms$Freq, nchar(all.terms$Var1))}
          all.terms<-all.terms[use.order,]
        
        ls.x<-lapply(ls.x, function(xx) xx[order(match(xx,all.terms$Var1))])
        return(sapply(ls.x, paste, collapse=' + '))
    }
#######        

stepAICc.df<-function(fits){
  fit.df<-rbind.fill(lapply(fits, function(x) as.data.frame(x[1:5])))
    null.deviance<-deviance(update(fits[[1]]$object, '.~1'))
  fit.df$r2.verhoef<-sapply(fits, function(x) {1-x$deviance/null.deviance})
  fit.df$r2.nakagawa<-sapply(fits, function(x){ 
      if(!is.null(x$object)){
        r2.nakagawa.glm(x$object)
      } else { r2.nakagawa.glm(update(fits[[1]]$object, x$formula))}})
  fit.df$k<-sapply(strsplit(paste(fit.df$formula),split='+',fixed=T), function(x){
                        x<-gsub('.*~','',x)
                        x<-gsub(' ','',x)
                        if(all(x!='1')) x<-c(x,'1')
                        length(x)})      
  fit.df$formula<-clean.formula(fit.df$formula)
  fit.df$dAIC<-fit.df$AIC-fit.df$AIC[1]
  fit.df$ModelID<-1:nrow(fit.df)
  fit.df<-subset(fit.df,select=c(deviance, df.resid,k,AIC,dAIC,r2.verhoef,r2.nakagawa,formula,ModelID))
  return(fit.df)}

#######

#######        

build.coef.df<-function(fitset){
  if(class(fitset[[1]])[1]=='list'){
        fitset<-lapply(fitset, function(x) x$object)
          fitset<-fitset[!sapply(fitset, is.null)]
  }

     if(class(fitset[[1]])[1]=='mer'){
        all.terms<-unique(unlist(lapply(fitset, function(x) names(x@fixef))))
         all.term.df<-rep(NA, length(all.terms))
         names(all.term.df)<-all.terms

          coef.df<-lapply(fitset, function(x) {
            coef.df<-all.term.df
            coef.df[match(names(x@fixef), all.terms,0)]=x@fixef
            return(coef.df)} )
     } else {      
      all.terms<-unique(unlist(lapply(fitset, function(x) names(coef(x)))))
        all.term.df<-rep(NA, length(all.terms))
        names(all.term.df)<-all.terms
    
      coef.df<-lapply(fitset, function(x) {
                coef.df<-all.term.df
                coef.df[match(names(coef(x)), all.terms,0)]=coef(x)
                return(coef.df)} )
      }
  as.data.frame(do.call('rbind',coef.df))
  }

    

          
  #######     