use "{file_link}\2019-1-30-024_FAHIM.dta", clear


gen lavgsal = log(avgsal)
gen lenrol = log(enrol)
gen lstaff = log(staff)

reg lavgsal bs lenrol lstaff lunch

** Assingment 1: Graphical Method **

quietly reg lavgsal bs lenrol lstaff lunch
predict uhat, resid

// plotting fitted value against residuals
rvfplot, yline(0)

// plotting individual independent variable against residuals
gen uhat_sq = uhat^2
scatter (uhat_sq bs), yline(0) name(bs)
scatter (uhat_sq lenrol), yline(0) name(enrol)
scatter (uhat_sq lstaff), yline(0) name(staff)
scatter (uhat_sq lunch), yline(0) name(lunch)
gr combine bs enrol staff lunch

** Assingment 2: Park test **

quietly reg lavgsal bs lenrol lstaff lunch
predict uhat, residual

gen uhat_sq = uhat^2
gen luhat_sq = log(uhat_sq)
gen lbs = log(bs)
gen llunch = log(lunch)

reg luhat_sq lbs lenrol lstaff lunch

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(4,.05)
scalar list  N R2 NR2 chi2


** Assingment 3: Glesjar test **

quietly reg lavgsal bs lenrol lstaff lunch
predict uhat, residual

gen abs_uhat = abs(uhat)

reg abs_uhat bs lenrol lstaff lunch

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(4,.05)
scalar list N R2 NR2 chi2

** Assingment 4: Gold field Quandt test **


// to select sub-samples

scalar c = 1848/6
scalar s_size = (1848 - c)/2
scalar list c s_size

sort bs

reg lavgsal bs lenrol lstaff lunch in 1/770
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

reg lavgsal bs lenrol lstaff lunch in 1079/1848
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

scalar FGQ = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,FGQ)
scalar list FGQ pvalue Fcrit

** Assingment 5: Breush-Pagan Godfrey test **
 
quietly reg lavgsal bs lenrol lstaff lunch
predict uhat, residual

gen uhat_sq = uhat^2
qui summarize uhat_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(1847-4)
scalar list RSS sigma_sq

gen p_i = uhat_sq/sigma_sq

reg p_i bs lenrol lstaff lunch

scalar ESS = 304.545126
scalar THETA = 0.5*ESS
scalar list THETA

scalar chi2=invchi2tail(4,.05)
scalar list THETA chi2

// STATA version Breush-Pagan Godfrey test
quietly reg lavgsal bs lenrol lstaff lunch
estat hettest

** Assingment 6: General White test **

quietly regress lavgsal bs lenrol lstaff lunch
predict uhat, residual

gen uhat_sq = uhat^2
gen bs_sq = bs^2
gen lenrol_sq = lenrol^2
gen lstaff_sq = lstaff^2
gen lunch_sq = lunch^2

reg uhat_sq bs lenrol lstaff lunch bs_sq lenrol_sq lstaff_sq lunch_sq

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(8,.05)
scalar list N R2 NR2 chi2

// STATA vestion white test
quietly regress lavgsal bs lenrol lstaff lunch
estat imtest, white

*** Autocorrelation ***

** Assingment 7: graphical method

quietly regress lavgsal bs lenrol lstaff lunch

gen time = _n
tsset time

predict uhat, residual
gen uhat_lag1 = L1.uhat
tsline uhat
twoway (scatter uhat_lag1 uhat)

** Assingment 8: Runs test **

quietly reg lavgsal bs lenrol lstaff lunch
predict uhat, residual

scatter uhat time, connect(l) yline(0) // to see the ups and down

runtest uhat, threshold(0)


** Assingment 9: Durbin Watson test

quietly regress lavgsal bs lenrol lstaff lunch
estat dwatson

** Assingment 10: Breusch-Godfrey test **

quietly regress lavgsal bs lenrol lstaff lunch
predict uhat, residuals

gen uhat_lag1 = L1.uhat
gen uhat_lag2 = L2.uhat

reg uhat uhat_lag1 uhat_lag2, noconstant

scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(2,.05)
scalar list N R2 NR2 chi2

// STATA version of Breusch-Godfrey test:

quietly regress lavgsal bs lenrol lstaff lunch
estat bgodfrey, lags(1/2)

*** Multicollinearity test ***

** Assingment 11:  R-squared and t value **

regress lavgsal bs lenrol lstaff lunch


* Pair-wise correlations

pwcorr bs lenrol lstaff lunch

** Assingment 12: Auxiliary regression for multicollinearity **

// bs
reg bs lenrol lstaff lunch

// lenrol
reg lenrol bs lstaff lunch

// lunch
reg lunch  lenrol bs lstaff 

// lstaff
reg lstaff lunch  lenrol bs

** Assingment 13: partial correlations **

pcorr lenrol lstaff lunch bs
pcorr lstaff bs lenrol lunch
pcorr lunch bs lenrol lstaff
pcorr bs lenrol lstaff lunch

** Assingment 14: condition index **

net install collin // collin is a third party library
collin bs lenrol lstaff lunch

** Assingment 15: Ramsey reset test

// Stata's version of the ramsey reset test
quietly regress lavgsal bs lenrol lstaff lunch
estat ovtest

// equivlent Mannual Ramsey Reset test

reg lavgsal bs lenrol lstaff lunch

predict yhat, xb
gen yhat_2 = yhat^2
gen yhat_3 = yhat^3
gen yhat_4 = yhat^4

reg lavgsal bs lenrol lstaff lunch yhat_2 yhat_3 yhat_4

scalar N=_result(1)
scalar R2_old= 0.4882
scalar R2_new = 0.4885 
scalar new_reg = 3
scalar no_param = 8
scalar NR2=N*R2

scalar F_value = ((R2_new - R2_old) / new_reg) / ((1 - R2_new) / (N - no_param))
scalar F_crit = invFtail(6,4,.05)
scalar P_value = Ftail(6,4,F_value)

scalar list R2_old R2_new N new_reg no_param F_value F_crit P_value