logit_HK <-function(y, X, pid){

# FITS CML Dynamic Logit (Honore Kyriazidou, 2000) 

# y       : vector of binary response variables of dimension nT x 1
#           (n=units, T=times)
# X       : a continuous covariate of dimension nT x 1 



    r = length(pid)
    label = unique(pid)
    n = length(label)
    X = as.matrix(X)
    k = ncol(X)
    eta = rep(0,k+1)
    TT = r/n
    pid = as.matrix(pid)

    sc = 0
    sigma = 8*n^(-0.2)
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
            x_i = X[pid==il]
            for(t in 2:(TT-2)){
                for(s in (t+1):(TT - 1)){
                    hk = ((y_i[t] + y_i[s]) == 1)
                    if(hk){
			
                      
                        v_i = c((x_i[t] - x_i[s]),(y_i[t-1] - y_i[s+1] + (y_i[t+1] - y_i[s-1])*(as.numeric(s-t>1))))
						
			ndx = v_i%*%eta
                        w = dnorm((x_i[t+1] - x_i[s+1])/sigma)
                        
                        lk = lk + w*(y_i[t]*ndx - log(1+exp(ndx)))
                        
                        d1 = w*(y_i[t] - exp(ndx)/(1+exp(ndx)))

                        if(is.nan(d1)){
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

        sc = colSums(scv)
        iJ = solve(J)
        
        eta = eta - iJ%*%sc
        
    }
    
    out = list(eta=eta,scv=scv,J=J)
    
}
