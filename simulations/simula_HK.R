## "Advances in maximum likelihood estimation of
##  fixed-effects binary panel data models"
## Francesco Valentini, Claudia Pigini, and Francesco Bartolucci
## SIMULATION STUDY - Results reported in Tables 1,2,3,5,6,7,8,9,10,11,12.

## # -------- ESTIMATORS -------------

## # INF : Infeasible (HK 2000)
## # ML  : Profile Likelihood
## # HK  : Honore & Kyriazidou, 2000
## # DJ  : Half-Panel Jackknife (Dahene and Jockmans, 2015)
## # QE  : Pseudo-Conditional (Bartolucci & Nigro, 2012)
## # MPL : Modified Profile Likelihood (Bartolucci et al., 2018)

require(panelMPL)
require(cquad)
source("HK-1.R")
source("dj.R")

# ---------------------------- SETUP -------------------------------- #
nit   = 1000
ssize = c(250,500,1000) # n. of subjects
time  = c(4,5,7,9,13)     # n. of occasions
ncov  = 1          # n. of covariates

gamma = c(0,0.25,0.5,1,2) # state dependence parameter
be = 1 # regression coefficient



# ------------------------ COVARIATES MODEL ------------------------- #

# Individual effects

st = pi/sqrt(3)                 # st. dev. x

# -------------------------- ESTIMATORS ---------------------------- #

INF = matrix(0,nit,ncov+1)
ML = matrix(0,nit,ncov+1) 
QE  = matrix(0,nit,ncov+1)
MPL = matrix(0,nit,ncov+1)
HK  = matrix(0,nit,ncov+1)
C = matrix(0,nit,ncov+1)
DJ = matrix(0,nit,ncov+1)

for(nu in ssize){
    for(TT in time){
        for(ga in gamma){
            filename = sprintf("N%g_T%g_ga%g_be%g_IT%g",nu ,TT, ga, be, nit)

            id = (1:nu)%x%rep(1,TT)
            idmpl = (1:nu)%x%rep(1,(TT-1))

            it = 0
            while(it < 1000){
               
                            # -------- GENERATE ARTIFICIAL DATA -------------- #

                x = xmpl = alpha = NULL
                s = 0
                for(i in 1:nu){
                    s = s+1
                    xi = as.vector(rnorm(TT, mean=0, sd = st))
                    x = c(x, xi)
                    xmpl = c(xmpl ,xi[2:TT])
                    alpha = c(alpha, mean(xi[1:4]))
                }
                X = as.matrix(x)
                Xmpl = as.matrix(xmpl)
                eta = c(be,ga)
                data = sim_panel_logit(id,alpha,X,eta,dyn=T)
                y = data$yv

                j = 0; yl = rep(0, nu*TT); y0 = rep(0,nu); ympll = ympl = NULL
                for(i in 1:nu){
                    j = j + 1
                    yl[j] = NA; y0[i] = y[j]                    
                    for(t in 2:TT){
                        j = j + 1
                        yl[j] = y[j-1]
                        ympl =  c(ympl,y[j])
                        ympll = c(ympll, y[j-1])
                    }
                }
                
                
                
                ALPHA = alpha%x%rep(1,(TT-1))
                
                # ------------ MODEL ESTIMATION ------------- #

                # HK
                out_HK = logit_HK(y,X,id)
                cond = out_HK$eta

                if(all(is.na(cond))){
                    print(0)
                }
                else{
                    
                    it = it + 1
                    HK[it,] = out_HK$eta


                    
                   
                                        # ML
                    out_mpl = panelMPL(ympl~Xmpl, panel=idmpl, model="dynLogit", y0=y0, method="BFGS", R = 200)
                    ML[it,] = out_mpl$mle
                    
                                        # MPL
                    MPL[it,] = out_mpl$est
                    
                    
                                        # QE
                    out_qe = cquad_pseudo(id,y,X)
                    QE[it,] = out_qe$coefficients
                    
                    
                                        # INFEASIBLE
                    out_inf = glm(ympl ~ -1 + Xmpl + ympll + ALPHA, family=binomial())
                    INF[it,] = out_inf$coefficients[1:2]

                    if(TT>6){
                                        #DJ
                    out_dj = dj(ympl,Xmpl,idmpl,y0)
                    DJ[it,] = out_dj$bh
                    }
 
                }


                        # --------------- OUTPUT ------------------- #

                if(it>1){
                    print("-----------------------------")
                    cat(c("it = ", it, "\n"))
                    cat(c("nu = ", nu, "\n"))
                    cat(c("TT = ", TT, "\n"))
                    cat(c("ga = ", ga, "\n"))
                    
                    true = matrix(c(1,ga),1,2)
                    cols = c("ML", "QE", "MPL", "HK", "INF","DJ")
                    rows = c("b1","ga")
                    MEAN = matrix(0,ncov+1,6)
                    MEAN[,1] = colMeans(ML[1:it,])
                    MEAN[,2] = colMeans(QE[1:it,])
                    MEAN[,3] = colMeans(MPL[1:it,])
                    MEAN[,4] = colMeans(HK[1:it,])
                    MEAN[,5] = colMeans(INF[1:it,])
                    MEAN[,6] = colMeans(DJ[1:it,])
                    colnames(MEAN) = cols
                    rownames(MEAN) = rows
                    
                    BIAS = matrix(0,ncov+1,6)
                    BIAS[,1] = colMeans(ML[1:it,]) - true
                    BIAS[,2] = colMeans(QE[1:it,]) - true
                    BIAS[,3] = colMeans(MPL[1:it,]) - true
                    BIAS[,4] = colMeans(HK[1:it,]) - true
                    BIAS[,5] = colMeans(INF[1:it,]) - true
                    BIAS[,6] = colMeans(DJ[1:it,]) - true
                    colnames(BIAS) = cols
                    rownames(BIAS) = rows
                    
                    true = matrix(rep(true,it),it,byrow=T)
                    MAE = matrix(0,ncov+1,6)
                    MAE[,1] = colMeans(abs(ML[1:it,] -true)) 
                    MAE[,2] = colMeans(abs(QE[1:it,] -true)) 
                    MAE[,3] = colMeans(abs(MPL[1:it,] -true))
                    MAE[,4] = colMeans(abs(HK[1:it,] -true)) 
                    MAE[,5] = colMeans(abs(INF[1:it,] -true))
                    MAE[,6] = colMeans(abs(DJ[1:it,] -true))

                    colnames(MAE) = cols
                    rownames(MAE) = rows
                    
                    RMSE = matrix(0,ncov+1,6)
                    RMSE[,1] = sqrt(colMeans((ML[1:it,] -true)^2))
                    RMSE[,2] = sqrt(colMeans((QE[1:it,] -true)^2))
                    RMSE[,3] = sqrt(colMeans((MPL[1:it,] -true)^2))
                    RMSE[,4] = sqrt(colMeans((HK[1:it,] -true)^2))
                    RMSE[,5] = sqrt(colMeans((INF[1:it,] -true)^2))
                    RMSE[,6] = sqrt(colMeans((DJ[1:it,] -true)^2))
                    
                    colnames(RMSE) = cols
                    rownames(RMSE) = rows

                    cat("Mean of regression coefficients\n")
                    print(MEAN)
                    cat("Mean bias\n")
                    print(BIAS)
                    cat("Mean absolute error\n")
                    print(MAE)
                    cat("Root mean square error\n")
                    print(RMSE)
                }
                       
                            
            }
            save.image(filename) 
            
        }
    }           
}

