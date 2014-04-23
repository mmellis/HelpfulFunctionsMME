# General subfunctions #
tryN<-function(f) tryCatch(f, error=function(e) return(NA))

unlogit<-function(x) {exp(x)/(1+exp(x))}

std<-function(x){ (x-mean(x,na.rm=T))/sd(x,na.rm=T)}

file.clean<-function(fname) {
  fn<-dir(pattern=fname, all.files=T, full.names=T)
  fn<-fn[!(grepl('rnw', fn)|grepl('tex',fn)|grepl('pdf',fn))]
  file.remove(fn)
  return(fn)
  }
  
# R^2 functions #  
r2.corr.mer <- function(m) {
   lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
   summary(lmfit)$r.squared
}
r2.corr.merPOP<-function(m) {
  lmfit<-lm(model.response(model.frame(m)) ~ (model.matrix(m) %*% lme4::fixef(m)))
  summary(lmfit)$r.squared
  }
create.null.model<-function(fit){
  form<-paste(formula(fit))
    if(class(fit)[1]=='mer') form[3]<-gsub('.*\\+ \\(','1 + (',form[3]) else form[3]<-1
  form<-formula(paste(form[2],form[3], sep='~'))  
    fit<-update(fit, form)
    return(fit)}
r2.verhoef<-function(m) {
  null.model<-create.null.model(m)
  1-deviance(m)/deviance(null.model)}  
  
  
r2.nakagawa.glm<-function(m, varType=2){
   if('glm' %in% class(m)){
        VarF<-var(as.vector(model.matrix(m) %*% coef(m)))
        VarF/(VarF + pi^2/3)
   } else if('glmerMod' %in% class(m)) {
        VarF<-var(as.vector(model.matrix(m) %*% lme4::fixef(m)))
        VarR<-c(det(lme4::VarCorr(m)[[1]]), sum(diag(lme4::VarCorr(m)[[1]])))
        list(marginal=c(VarF/(VarF+VarR+pi^2/3), VarF/(var(fitted(m))+pi^2/3))[varType],
             conditional=c((VarF+VarR)/(VarF+VarR+pi^2/3),var(fitted(m))/(var(fitted(m))+pi^2/3))[varType]) 
   }else { return(NA) }
   }
   
rename.var<-function(dataframe, from.name, to.name){
   names(dataframe)[match(from.name, names(dataframe))]=to.name    
   return(dataframe)}

last<-function(x,na.rm=T){
  if(na.rm) x<-x[!is.na(x)]
  return(rev(x)[1])
}

substr1<-function(x,ch=1){  
  unlist(lapply(strsplit(paste(x),split=''), 
    function(xx) paste(xx[ch],collapse=''))) }
    
fill.column<-function(x){
  values<-which(!is.na(x))
  rep(x[values],times=diff(c(values,length(x)+1)))}
  
