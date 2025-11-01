gof_bin <- function(B, bimodel = NULL, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if( is.null(bimodel)) stop("model object is not given.")
  y <- bimodel$model[, 1]
  n <- length(y)
  p1f <- bimodel$fitted.values
  dfx <- model.matrix(bimodel)
  fam <- bimodel$family
  como <- combn(n,2)
  ind1 <- como[1,]
  ind2 <- como[2,]


  ## F(Y), F(Y-1), P(Y)

  Fy <- ifelse(y==1,1,1-p1f)
  Fy1 <-  ifelse(y==0,0,1-p1f)
  dy <- ifelse(y==1,p1f,1-p1f)

  ## The test statistic is calculated by parts

  int1new <- sum(1/3*(Fy^4+Fy1^4-Fy^3*Fy1-Fy*Fy1^3-2*Fy^3-Fy1^3+3*Fy^2*Fy1+Fy^2-2*Fy*Fy1+Fy1^2)/dy^2)


  Fy1new <- sort(Fy1)
  Fynew <- Fy[sort(Fy1,index=TRUE)$ix]
  dynew <- dy[sort(Fy1,index=TRUE)$ix]




  ###no overlap
  group1 <- which(Fy1new[ind2]>Fynew[ind1])
  ## overlap, and keep order. k=ind1,k'=ind2,m=ind2,m'=ind1
  group2 <- which(Fy1new[ind2]<=Fynew[ind1]&Fynew[ind2]>Fynew[ind1])
  ## overlap, and switch order for on right side. k=ind1,k'=ind2,m=ind1,m'=ind2
  group3 <- which(Fy1new[ind2]<=Fynew[ind1]&Fynew[ind2]<=Fynew[ind1])
  group4 <- which(Fy1new[ind2]<=Fynew[ind1])

  ind2new1 <- sum(((Fynew^3-Fy1new^3)/6/dynew)[ind1[group1]])

  ind2new2 <- sum(((Fynew^3-Fy1new^3-3*Fynew^2+3*Fy1new^2+2*Fynew-2*Fy1new)/dynew/6)[ind2[group1]])


  indnew3p1 <- sum((-Fy1new[ind1[group4]]^3-2*Fy1new[ind2[group4]]^3*(1-dynew[ind1[group4]])+3*Fy1new[ind1[group4]]*Fy1new[ind2[group4]]^2)/6/dynew[ind1[group4]])



  indnew3p2 <- sum((1-dynew[ind1[group2]])*(1-dynew[ind2[group2]])/dynew[ind1[group2]]/dynew[ind2[group2]]*(Fynew[ind1[group2]]^3-Fy1new[ind2[group2]]^3)/3) +
    sum((1-dynew[ind1[group3]])*(1-dynew[ind2[group3]])/dynew[ind1[group3]]/dynew[ind2[group3]]*(Fynew[ind2[group3]]^3-Fy1new[ind2[group3]]^3)/3) -
    sum(((1-dynew[ind1[group2]])*(Fy1new[ind2[group2]])+(1-dynew[ind2[group2]])*(Fy1new[ind1[group2]]))/dynew[ind1[group2]]/dynew[ind2[group2]]*(Fynew[ind1[group2]]^2-Fy1new[ind2[group2]]^2)/2) -
    sum(((1-dynew[ind1[group3]])*(Fy1new[ind2[group3]])+(1-dynew[ind2[group3]])*(Fy1new[ind1[group3]]))/dynew[ind1[group3]]/dynew[ind2[group3]]*(Fynew[ind2[group3]]^2-Fy1new[ind2[group3]]^2)/2) +
    sum((Fy1new[ind1[group2]])*(Fy1new[ind2[group2]])/dynew[ind1[group2]]/dynew[ind2[group2]]*(Fynew[ind1[group2]]-Fy1new[ind2[group2]])) +
    sum((Fy1new[ind1[group3]])*(Fy1new[ind2[group3]])/dynew[ind1[group3]]/dynew[ind2[group3]]*(Fynew[ind2[group3]]-Fy1new[ind2[group3]]))

  indnew3p3 <- sum(1/6/dynew[ind2[group2]]*(Fynew[ind2[group2]]^3+2*Fynew[ind1[group2]]^3-2*Fynew[ind1[group2]]^3*Fynew[ind2[group2]]+
                                              2*Fynew[ind1[group2]]^3*Fy1new[ind2[group2]]-
                                              3*Fynew[ind1[group2]]^2-3*Fynew[ind2[group2]]^2+3*Fynew[ind1[group2]]^2*Fynew[ind2[group2]]-
                                              6*Fynew[ind1[group2]]^2*Fy1new[ind2[group2]]+6*Fynew[ind1[group2]]*Fy1new[ind2[group2]]+
                                              2*Fynew[ind2[group2]]-2*Fy1new[ind2[group2]]))+
    sum(1/6/dynew[ind1[group3]]*(Fynew[ind1[group3]]^3+2*Fynew[ind2[group3]]^3-2*Fynew[ind2[group3]]^3*Fynew[ind1[group3]]+
                                   2*Fynew[ind2[group3]]^3*Fy1new[ind1[group3]]-
                                   3*Fynew[ind2[group3]]^2-3*Fynew[ind1[group3]]^2+3*Fynew[ind2[group3]]^2*Fynew[ind1[group3]]-
                                   6*Fynew[ind2[group3]]^2*Fy1new[ind1[group3]]+6*Fynew[ind2[group3]]*Fy1new[ind1[group3]]+
                                   2*Fynew[ind1[group3]]-2*Fy1new[ind1[group3]]))



  ## Test statistic
  disnull <- int1new+2*(ind2new1+ind2new2+indnew3p1+indnew3p2+indnew3p3)

  ## Bootstrap
  disin <- rep(0,B)

  for(repin in 1:B){

    yr <- rbinom(n,size=1,p1f)
    #####!!!changes made
    #####!!!
    dfr <- data.frame(yr=yr,dfx[,-1])

    #fit marginal model
    bimodelr <- glm(yr~.,family=fam,data=dfr)
    p1fr <- bimodelr$fitted.values
    #####!!!
    #####!!!




    Fyr <- ifelse(yr==1,1,1-p1fr)
    Fy1r <-  ifelse(yr==0,0,1-p1fr)
    dyr <- ifelse(yr==1,p1fr,1-p1fr)

    int1newr <- sum(1/3*(Fyr^4+Fy1r^4-Fyr^3*Fy1r-Fyr*Fy1r^3-2*Fyr^3-Fy1r^3+3*Fyr^2*Fy1r+Fyr^2-2*Fyr*Fy1r+Fy1r^2)/dyr^2)



    Fy1new <- sort(Fy1r)
    Fynew <- Fyr[sort(Fy1r,index=TRUE)$ix]
    dynew <- dyr[sort(Fy1r,index=TRUE)$ix]



    group1 <- which(Fy1new[ind2]>Fynew[ind1])
    ## overlap, and keep order. k=ind1,k'=ind2,m=ind2,m'=ind1
    group2 <- which(Fy1new[ind2]<=Fynew[ind1]&Fynew[ind2]>Fynew[ind1])
    ## overlap, and switch order for on right side. k=ind1,k'=ind2,m=ind1,m'=ind2
    group3 <- which(Fy1new[ind2]<=Fynew[ind1]&Fynew[ind2]<=Fynew[ind1])
    group4 <- which(Fy1new[ind2]<=Fynew[ind1])

    ind2new1r <- sum(((Fynew^3-Fy1new^3)/6/dynew)[ind1[group1]])

    ind2new2r <- sum(((Fynew^3-Fy1new^3-3*Fynew^2+3*Fy1new^2+2*Fynew-2*Fy1new)/dynew/6)[ind2[group1]])


    indnew3p1r <- sum((-Fy1new[ind1[group4]]^3-2*Fy1new[ind2[group4]]^3*(1-dynew[ind1[group4]])+3*Fy1new[ind1[group4]]*Fy1new[ind2[group4]]^2)/6/dynew[ind1[group4]])



    indnew3p2r <- sum((1-dynew[ind1[group2]])*(1-dynew[ind2[group2]])/dynew[ind1[group2]]/dynew[ind2[group2]]*(Fynew[ind1[group2]]^3-Fy1new[ind2[group2]]^3)/3) +
      sum((1-dynew[ind1[group3]])*(1-dynew[ind2[group3]])/dynew[ind1[group3]]/dynew[ind2[group3]]*(Fynew[ind2[group3]]^3-Fy1new[ind2[group3]]^3)/3) -
      sum(((1-dynew[ind1[group2]])*(Fy1new[ind2[group2]])+(1-dynew[ind2[group2]])*(Fy1new[ind1[group2]]))/dynew[ind1[group2]]/dynew[ind2[group2]]*(Fynew[ind1[group2]]^2-Fy1new[ind2[group2]]^2)/2) -
      sum(((1-dynew[ind1[group3]])*(Fy1new[ind2[group3]])+(1-dynew[ind2[group3]])*(Fy1new[ind1[group3]]))/dynew[ind1[group3]]/dynew[ind2[group3]]*(Fynew[ind2[group3]]^2-Fy1new[ind2[group3]]^2)/2) +
      sum((Fy1new[ind1[group2]])*(Fy1new[ind2[group2]])/dynew[ind1[group2]]/dynew[ind2[group2]]*(Fynew[ind1[group2]]-Fy1new[ind2[group2]])) +
      sum((Fy1new[ind1[group3]])*(Fy1new[ind2[group3]])/dynew[ind1[group3]]/dynew[ind2[group3]]*(Fynew[ind2[group3]]-Fy1new[ind2[group3]]))

    indnew3p3r <- sum(1/6/dynew[ind2[group2]]*(Fynew[ind2[group2]]^3+2*Fynew[ind1[group2]]^3-2*Fynew[ind1[group2]]^3*Fynew[ind2[group2]]+
                                                 2*Fynew[ind1[group2]]^3*Fy1new[ind2[group2]]-
                                                 3*Fynew[ind1[group2]]^2-3*Fynew[ind2[group2]]^2+3*Fynew[ind1[group2]]^2*Fynew[ind2[group2]]-
                                                 6*Fynew[ind1[group2]]^2*Fy1new[ind2[group2]]+6*Fynew[ind1[group2]]*Fy1new[ind2[group2]]+
                                                 2*Fynew[ind2[group2]]-2*Fy1new[ind2[group2]]))+
      sum(1/6/dynew[ind1[group3]]*(Fynew[ind1[group3]]^3+2*Fynew[ind2[group3]]^3-2*Fynew[ind2[group3]]^3*Fynew[ind1[group3]]+
                                     2*Fynew[ind2[group3]]^3*Fy1new[ind1[group3]]-
                                     3*Fynew[ind2[group3]]^2-3*Fynew[ind1[group3]]^2+3*Fynew[ind2[group3]]^2*Fynew[ind1[group3]]-
                                     6*Fynew[ind2[group3]]^2*Fy1new[ind1[group3]]+6*Fynew[ind2[group3]]*Fy1new[ind1[group3]]+
                                     2*Fynew[ind1[group3]]-2*Fy1new[ind1[group3]]))



    disin[repin] <-int1newr+2*(ind2new1r+ind2new2r+indnew3p1r+indnew3p2r+indnew3p3r)

  }
  pvalue <- length(which(disin>disnull))/B
  return(list(test_stat = disnull, p_value = pvalue))
}
