logit_HK <-function(y, Xc, Xd, pid,c=8,conv=(-1/5)){

    ## FITS CML Dynamic Logit (Honore Kyriazidou, 2000) 
    ## y  = Dependent Variable
    ## Xc = Continuous Covariates
    ## Xd = Discrete Covariates
    ## pid = Index
    
    r = length(pid)
    label = unique(pid)
    n = length(label)
    X = as.matrix(cbind(Xc,Xd))
    k = ncol(X)

    
    eta = rep(0,k+1)
                                     

    Tv = rep(0,n)
    for(i in 1:n){
        Tv[i] = length(y[pid==label[i]])
    }

    
    sc = 0
    sigma = c*n^(conv)*diag(k)
    it = 0; lk = 0; lk0 = -Inf
    while(any(abs(sc)>10^-6) | it==0){
        it = it+1
        lk0 = lk
        lk = 0
        J = 0
        scv = matrix(0,n,k+1)
        for(i in 1:n){
            il = label[i]
            y_i = y[pid==il] 
            x_i = as.matrix(X[pid==il,])
            xd_i = as.matrix(Xd[pid==il,])
            if(Tv[i]>3){
                for(t in 2:(Tv[i]-2)){
                    for(s in (t+1):(Tv[i] - 1)){
                        
                        hk = (all((xd_i[t+1,] - xd_i[s+1,]) ==0 ))*((y_i[t] + y_i[s]) == 1)
                        
                        if(hk){
                            
                            v_i = c((x_i[t,1:k] - x_i[s,1:k]),(y_i[t-1] - y_i[s+1] + (y_i[t+1] - y_i[s-1])*(as.numeric(s-t>1))))
                            
                            ndx = v_i%*%eta
                            w = dmvnorm((x_i[t+1,1:k] - x_i[s+1,1:k])%*%solve(sigma))
                            
                                
                            lk = lk + w*(y_i[t]*ndx - log(1+exp(ndx)))
                            
                            d1 = w*(y_i[t] - exp(ndx)/(1+exp(ndx)))
                            
                            if(any(is.nan(d1))){
                                out = list(eta=NA,scv=NA,J=NA)
                                return(out)
                                
                            }                       
                            else{
                                
                                
                                scv[i,] =  scv[i,] + d1%*%t(v_i)
                                
                                d2 = as.vector(w*(exp(ndx)/(1+exp(ndx))^2))
                                J = J - d2*(v_i%o%v_i)
                            }
                            
                        }
                        
                        
                    }
                }
            }
        }
        sc = colSums(scv)
        iJ = solve(J)
        
        eta = eta - iJ%*%sc
    }

    out = list(eta=eta,scv=scv,J=J)
    
}

