dyn_logit_ml  <-  function(id,y,X){

    label = unique(id); n = length(label) 
    X = as.matrix(X); k = ncol(X)
    Tv = rep(0, n)
    for (i in 1:n) Tv[i] = sum(id == label[i])
    TT = max(Tv)
    
    ## Prepare data
    Xal = diag(n)%x%rep(1,TT)
    indc = indr = NULL
    for(i in 1:n){
        y_i = y[id==label[i]]
        sui = sum(y_i)
        if(sui==0 | sui==TT){
            indc = c(indc, i)
            indr = c(indr, rep(1,Tv[i]))
        }else{
            indr = c(indr, rep(0,Tv[i]))
        }
    }
    Xal = Xal[!indr, ]
    if(!is.null(indc)) Xal = Xal[,-indc]
    NT = sum(Tv); NTq = length(id[!indr])

    
     ## MAXIMUM LIKELIHOOD ESTIMATION
    y1 = y[!indr]; X1 = X[!indr,]
    X1 = cbind(X1, Xal); npar = ncol(X1)
    data = as.data.frame(cbind(y1,X1))
    mod = glm(y1 ~ . -1, data = data, family = binomial(link="logit"))
    be = mod$coefficients[1:k]; al = mod$coefficients[(k+1):npar]    

    ## MLE variance
    AL = rep(NA,NT)
    scb = matrix(0,n,k)
    sca = matrix(0,n,ncol(Xal)); j = 0
    for(i in 1:n){
        il = label[i]
        y_i = y[id == il]; sui = sum(y_i)
        if(sui>0 & sui<Tv[i]){
            j = j+1
            x_i = as.matrix(X[id == il,])
            xb = x_i%*%be + al[j]
            F_i = as.vector(exp(xb)/(1 + exp(xb)))
            scb[i,] = t(y_i - F_i)%*%x_i
            sca[i,j] = sum(y_i - F_i)
            AL[id == i] = al[j]
        }
    }
    scv = cbind(scb,sca)
    J = summary(mod)$cov.scaled
    Vr = J%*%(t(scv)%*%scv)%*%J
    se_be = sqrt(diag(as.matrix(Vr[1:k,1:k])))


    out = list(be = be, se_be = se_be, Vr = Vr)
}
