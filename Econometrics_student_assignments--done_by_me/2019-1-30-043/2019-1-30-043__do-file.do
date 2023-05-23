use "link\2019-1-30-043.dta", clear

regress lwage looks educ exper female

*** Assignment 1: raphical method 

quietly regress lwage looks educ exper female
predict epsilon, resid

rvfplot, yline(0)

gen epsilon_sq = epsilon^2
scatter (epsilon_sq looks), yline(0) name(looks)
scatter (epsilon_sq educ), yline(0) name(educ)
scatter (epsilon_sq exper), yline(0) name(exper)
scatter (epsilon_sq female), yline(0) name(female)
graph combine looks educ exper female


*** Assignment 2: Park test

quietly regress lwage looks educ exper female
predict epsilon, residual
gen epsilon_sq = epsilon^2

gen lepsilon_sq = log(epsilon_sq)
gen llooks = log(looks)
gen lexper = log(exper)
gen leduc = log(educ)

regress lepsilon_sq llooks leduc lexper female

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(4,.05)
scalar list N R2 LM crit_val

*** Assignment 3: Glesjer test

quietly regress lwage looks educ exper female
predict epsilon, residual

gen abs_epsilon = abs(epsilon)

regress abs_epsilon looks educ exper female

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(4,.05)
scalar list N R2 LM crit_val

*** Assignment 4: Gold field Quandt test

*sub samples*

scalar c = 1260/5
scalar sample_size = (1260 - c)/2
scalar list c sample_size

sort exper

regress lwage looks educ exper female in 1/504
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

regress lwage looks educ exper female in 757/1260
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit


*** Assignment 5:  Breush-Pagan Godfrey test
 
quietly regress lwage looks educ exper female
predict epsilon, residual

gen epsilon_sq = epsilon^2
qui summarize epsilon_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(1260-4)
scalar list RSS sigma_sq

gen p_i = epsilon_sq/sigma_sq
regress p_i looks educ exper female

scalar ess = 6.08737094
scalar theta = 0.5*ess

scalar crit_val=invchi2tail(4,.05)
scalar list theta crit_val

//Stata version of BPG test
estat hettest

*** Assignment 6:  White test

quietly regress lwage looks educ exper female
predict epsilon, residual

gen epsilon_sq = epsilon^2
gen looks_sq = looks^2
gen educ_sq = educ^2
gen exper_sq = exper^2

regress epsilon_sq looks educ exper female looks_sq educ_sq exper_sq

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(7,.05)
scalar list N R2 LM crit_val

//STATA white test
estat imtest, wite

*** Autocorrelation test ***

*** Assignment 7:  graphical method

quietly regress lwage looks educ exper female
predict epsilon, residual

gen time = _n
tsset time

gen epsilon_lag1 = L1.epsilon
twoway (scatter epsilon_lag1 epsilon) (lfit epsilon_lag1 epsilon)


*** Assignment 8:  Runs test

quietly regress lwage looks educ exper female
predict epsilon, residual

scatter epsilon time, connect(l) yline(0)
runtest epsilon, threshold(0)


*** Assignment 9:  Durbin Watson test -- by stata

quietly regress lwage looks educ exper female
estat dwatson

*** Assignment 10:  Breusch-Godfrey test

quietly regress lwage looks educ exper female
estat bgodfrey, lags(1/2)

*or

quietly regress lwage looks educ exper female
predict epsilon, residuals

gen epsilon_lag1 = L1.epsilon
gen epsilon_lag2 = L2.epsilon

regress epsilon epsilon_lag1 epsilon_lag2, noconstant

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2_crit=invchi2tail(2,.05)
scalar list N R2 LM chi2_crit


***** Multicollinearity test ******

*** Assignment 11: looking at R-squared and t value

regress lwage looks educ exper female


*** Assignment 12: Pair-wise correlations

pwcorr looks educ exper female


*** Assignment 13:  Auxiliary regress for multicollinearity

// for looks
regress looks educ exper female

// for educ
regress educ looks exper female

// for female
regress female  educ looks exper 

// for exper
regress exper female  educ looks


*** Assignment 14: partial correlations

pcorr looks educ exper female
pcorr educ exper female looks
pcorr exper looks educ female
pcorr female looks educ exper

*** Assignment 15: condition index

net install collin // third parrty library
collin looks educ exper female


*** Assignment 16:  Ramsey reset test

// Stata's version of the ramsey reset test
quietly regress lwage looks educ exper female
ovtest


// equivlent Mannual Reset test
regress lwage looks educ exper female

predict yhat, xb
gen yhat2 = yhat^2
gen yhat3 = yhat^3
gen yhat4=yhat^4

regress lwage looks educ exper female yhat2 yhat3 yhat4

scalar N=_result(1)
scalar R2_old= 0.3360
scalar R2_new = 0.3373
scalar new_var = 3
scalar tot_param = 8

scalar F_value = ((R2_new - R2_old) / new_var) / ((1 - R2_new) / (N - tot_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)

scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value