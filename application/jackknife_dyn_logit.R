hpjl  <- function(id,y,X){

    # LAGGED Y must be in X
    
    r = length(id)
    label = unique(id)
    n = length(label)
    TT = r/n
    pt = rep(1,n)%x%(1:TT)
    labelt = unique(pt)
    k = ncol(as.matrix(X))
    
    
                                        # Half-panels
    id1 = id[labelt<=TT/2]
    id2 = id[labelt>TT/2]
    
    y1 = y[labelt<=TT/2]
    y2 = y[labelt>TT/2]

    X1 = X[labelt<=TT/2,]
    X2 = X[labelt>TT/2,]

    
    out_full = dyn_logit_ml(id,y,X)
    b = out_full$be 
    
    out_sub1 = dyn_logit_ml(id1,y1,X1)
    b1 = out_sub1$be #first sub-sample
    V1 = out_sub1$Vr
    
    
    out_sub2 = dyn_logit_ml(id2,y2,X2)
    b2 = out_sub2$be #second-subsample
    V2 = out_sub2$Vr

    
    bh = (2*b) - (b1 + b2)*0.5 #jackknife estimator
    VV = (V1[1:k,1:k] + V2[1:k,1:k])*0.5
    se_bh = sqrt(diag(VV))

    out = list(bh = bh,b = b,b1 = b1,b2 = b2, V1 = V1, V2 = V2, VH = VV,se_bh = se_bh )
}
