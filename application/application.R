## "Advances in maximum likelihood estimation of
##  fixed-effects binary panel data models"
## Francesco Valentini, Claudia Pigini, and Francesco Bartolucci
## APPLICATION - Results reported in Table 4


require(foreign)
require(panelMPL)
require(cquad)
source("HK-3.R")
source("jackknife_dyn_logit.R")
source("logit-ml.R") 
require(mvtnorm)

set.seed(12454845)

## DATA PSID 1979-1985 
psid = read.dta("psid_example.dta")

## FULL ID
id = as.matrix(psid$id2)

## EXPLANATORY VARIABLES
X = as.matrix(cbind(psid$age,psid$age2,psid$tempinc,psid$kd2,psid$kd5,psid$kd17))
y = psid$part

Xc = as.matrix(cbind(psid$age,psid$age2,psid$tempinc))
Xd = as.matrix(cbind(psid$kd2,psid$kd5,psid$kd17))



## EXCLUDE 1979 
XX = X[psid$year>79,] 

## DEPENDENT VARIABLE 
yy = psid$part[psid$year>79]

## INITIAL CONDITIONS FOR MPL 
y0 = psid$part[psid$year==79]

## LAG DEP VARIABLE
yl = psid$lag_part[psid$year>79]

## INCLUDE LAGGED PART. IN COVARIATES
XXX = cbind(XX,yl)

## RESTRICTED ID
idx = id[psid$year>79]



## PSEUDO CML
out_pseudo = cquad_pseudo(id,y,X)
pseudo_coeff = out_pseudo$coefficients

## HONORÈ & KYRIAZIDOU
out_HK = logit_HK(y,Xc,Xd,id,c=128)
hk_coeff = out_HK$eta
se_hk = sqrt(diag((solve(out_HK$J)%*%(t(out_HK$scv)%*%out_HK$scv)%*%solve(out_HK$J))))



## MAXIMUM LIKELIHOOD
out_ml = dyn_logit_ml(idx,yy,XXX)    
ml_coeff= out_ml$be
ml_se = out_ml$se_be

##MODIFIED PROFILE LIKELIHOOD
out_mpl = panelMPL(yy~XX, panel=idx, model="dynLogit", y0=y0, method="BFGS", R = 200)
mpl_coeff = out_mpl$est


## JACKKNIFE
out_dj = hpjl(idx,yy,XXX)
dj_coeff = out_dj$bh


## ######## PRINT ESTIMATION RESULTS ################

covariates = rbind("age","age2","tempincome","kd2","kd5","kd17","y_lag")
colms = rbind("coeff","s.e.","t","p-value")

print("ML")
ML = cbind(out_ml$be,out_ml$se_be,(out_ml$be/out_ml$se_be),2*pnorm(abs(out_ml$be/out_ml$se_be),lower.tail=FALSE))
rownames(ML) = covariates
colnames(ML) = colms
print(ML)

print("PSEUDO")
PSEUDO = cbind(out_pseudo$coefficients,out_pseudo$ser,(out_pseudo$coefficients/out_pseudo$ser),2*pnorm(abs(out_pseudo$coefficients/out_pseudo$ser),lower.tail=FALSE))
rownames(PSEUDO) = covariates
colnames(PSEUDO) = colms
print(PSEUDO)


print("HK")
HK = cbind(out_HK$eta,se_hk,(hk_coeff/se_hk),2*pnorm(abs((hk_coeff/se_hk)),lower.tail=FALSE))
rownames(HK) = covariates
colnames(HK) = colms
print(HK)


print("MPL")
MPL = cbind(out_mpl$est,sqrt(diag(solve(-out_mpl$hessian))),(out_mpl$est/sqrt(diag(solve(-out_mpl$hessian)))),2*pnorm(abs(out_mpl$est/sqrt(diag(solve(-out_mpl$hessian)))),lower.tail=FALSE))
rownames(MPL) = covariates
colnames(MPL) = colms
print(MPL)

print("JACKKNIFE")
JK = cbind(out_dj$bh,out_dj$se_bh,(out_dj$bh/out_dj$se_bh),2*pnorm(abs(out_dj$bh/out_dj$se_bh),lower.tail=FALSE))
rownames(JK) = covariates
colnames(JK) = colms
print(JK)

