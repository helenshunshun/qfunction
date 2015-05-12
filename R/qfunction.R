#' find out indicator for a2 for each individual
#' 
#' It is the best treatment choice at time 2, which maximizes Q2 function and get V2 function.
#' @param outcome is CD4 # at the end
#' @param tx1,tx2 are treatment choices at each decision point
#' @param cov1,cov2 are two covariate of each individual
#' @return indicator of a2, if value>0, a2=1,taking a2 treatment; if value<0, a2=0,not taking a2 
#' @export 
indicat<-function(outcome,tx1,tx2,cov1,cov2){ 
  beta<-summary(lm(outcome~tx1+tx2+cov1+cov2+tx2*cov1+tx2*cov2+tx1*tx2))
  psi<-round(beta[[4]][c(3,6,7,8)],2)
  value<-psi[2]*cov1+psi[3]*cov2+psi[4]*tx2+rep(psi[1],length(tx2))
  indi<-as.data.frame(cbind(value,indi= ifelse(value>0,1,0)))
  return(indi)
}

#' get maximum of outcome
#' 
#' @param outcome is CD4 # at the end
#' @param tx1,tx2 are treatment choices at each decision point
#' @param cov1,cov2 are two covariate of each individual
#' @return maximum value of outcome
#' @export
Q2<-function(outcome,tx1,tx2,cov1,cov2){
  gg<-indicat(outcome,tx1,tx2,cov1,cov2)
  tdata<-as.data.frame(cbind(outcome,tx1,tx2,cov1,cov2,gg[[2]]))
  m<-nrow(tdata)
  f<-round(gg[[1]][1:8],2)
  q2<-rep(0,m) 
  for (i in 1:m){
    if(tdata$indi[i]==1){
      q2[i]<-f[1]+f[4]*cov1[i]+f[5]*cov2[i]+f[2]*tx1[i]
      +f[3]*tx2[i]+f[6]*cov1[i]*tx2[i] +f[7]*cov2[i]*tx2[i]+f[8]*tx1[i]*tx2[i] 
      ## if >0, take A2 as 1
    } 
    if(tdata$indi[i]==0){
      q2[i]<-f[1]+f[4]*cov1[i]+f[5]*cov2[i]+f[2]*tx1[i] ## if <0, take A2 as 0 
    }
  }
  return (q2) 
}

#' find out indicator for a1 for each individual
#' 
#' It is the best treatment choice at time 1, which maximizes Q1 function and get V1 function.
#' @param outcome is CD4 # at the end
#' @param tx1,tx2 are treatment choices at each decision point
#' @param cov1,cov2 are two covariate of each individual
#' @return maximum value of outcome
#' @export
indicat.b<-function(outcome,tx1,tx2,cov1,cov2){ 
  v2<-Q2(outcome,tx1,tx2,cov1,cov2)
  beta.b<-summary(lm(v2~cov1+cov2+tx1+cov1*tx1+cov2*tx1))
  coefficient.b<-beta.b[[4]]
  psi.b<-round(beta.b[[4]][c(4,5,6)],2)
  value.b<-psi.b[2]*cov1+psi.b[3]*cov2+rep(psi.b[1],length(tx1))
  indi.b<-as.data.frame(cbind(value.b,indi.b= ifelse(value.b>0,1,0)))
  return(list(coefficient.b,indi.b))
}

#' get maximum of V2
#' 
#' This step is try to get back the trial history
#' @param outcome is CD4 # at the end
#' @param tx1,tx2 are treatment choices at each decision point
#' @param cov1,cov2 are two covariate of each individual
#' @return maximum value of outcome
#' @export
Q1<-function(outcome,tx1,tx2,cov1,cov2){
  gg.b<-indicat.b(outcome,tx1,tx2,cov1,cov2)
  tdata.b<-as.data.frame(cbind(outcome,tx1,tx2,cov1,cov2,gg.b[[2]]))
  m.b<-nrow(tdata.b)
  f.b<-round(gg.b[[1]][1:6],2)
  q1<-rep(0,m.b)
  for (i in 1:m.b){
    if(tdata.b$indi.b[i]==1){ 
      q1[i]<-f.b[1]+f.b[2]*cov1[i]+f.b[3]*cov2[i]+f.b[4]*tx1[i]+
        f.b[5]*cov1[i]*tx1[i]+f.b[6]*cov2[i]*tx1[i] ## if >0, take A1 as 1
    } 
    if(tdata.b$indi.b[i]==0){
      q1[i]<-f.b[1]+f.b[2]*cov1[i]+f.b[3]*cov2[i]+f.b[4]*tx1[i] ## if <0, take A1 as 0
    }
  }
  return (q1)
}

