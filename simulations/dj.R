# HALF PANEL JACKKNIFE ESTIMATOR
# DAHENE AND JOCKMANS, 2015

dj <- function(y, X, pid, y0){

    
    r = length(pid)
    label = unique(pid)
    n = length(label)
    TT = r/n
    pt = rep(1,n)%x%(1:TT)
    labelt = unique(pt)
    
    
                                        #id for half-panels
    hpid  = (1:n) %x% rep(1,(TT/2))
    
                                        # Half-panels
    
    
    y1 = y[labelt<=TT/2]
    y2 = y[labelt>TT/2]
    y20 = y[labelt==(TT/2)]
    X1 = X[labelt<=TT/2,]
    X2 = X[labelt>TT/2,]

    out_full = panelMPL(y~X, panel=pid, model="dynLogit", y0=y0, method="BFGS", R=10)
    b = out_full$mle
    
    out_sub1 = panelMPL(y1~X1, panel=hpid, model="dynLogit", y0=y0, method="BFGS", R=10)
    b1 = out_sub1$mle #first sub-sample

    
    
    out_sub2 = panelMPL(y2~X2, panel=hpid, model="dynLogit", y0=y20, method="BFGS", R=10)
    b2 = out_sub2$mle #second-subsample

    

    bh = 2*b - (b1 + b2)/2 #jackknife estimator


    

   
    out = list(bh = bh,b = b,b1 = b1,b2 = b2)
    
}


    
