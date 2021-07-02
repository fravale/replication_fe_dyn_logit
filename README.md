# replication_fe_dyn_logit

"Advances in maximum likelihood estimation of fixed-effects binary panel data models" 
Francesco Valentini, Claudia Pigini, and Francesco Bartolucci, (forthcoming).

Supplemtary material: the R code is split in two folders.

1) simulations/ simula_HK.R <-- reproduces the simulation design Honoré and Kyriazidou (2000) and the relative estimation results (Tabs 1-3,5-12);
                HK-1.R      <-- performs the Honoré and Kyriazidou (2000) estimator for the DL model with 1 continuous explanatory variable;
		            dj.R        <-- performs the half-panl jackknife proposed by Dhaene and Jochmans (2015).

2) application/ application.R         <-- reproduces estimation results reported in Table 4;
                logit-ml.R            <-- performs the ML estimator for the logit model with fixed-effects;
		            jackknife_dyn_logit.R <-- performs the half-panl jackknife proposed by Dhaene and Jochmans (2015) (based on logit-ml.R);
		            HK-3.R                <-- performs the Honoré and Kyriazidou (2000) estimator for the DL model with continuous and discrete explanatory variables;
	 	            psid_example.dta      <-- the dataset

Dependencies / additional R packages: cquad    <-- https://cran.r-project.org/package=cquad
	       		    	                    foreign  <-- https://cran.r-project.org/package=foreign
				                              mvtnorm  <-- https://cran.r-project.org/package=mvtnorm
                                      panelMPL <-- https://ruggerobellio.weebly.com/uploads/5/1/5/0/51505127/panelmpl_0.23.tar.gz
