use "C:\Users\nafew\OneDrive - East West University\Academic\1. Fall-2022\A-465\Oboni\2019-3-30-029.DTA", clear



reg educ motheduc fatheduc abil urban

// graphical method 

quietly reg educ motheduc fatheduc abil urban
predict uhat, resid
rvfplot, yline(0)

gen uhat_sq = uhat^2
scatter (uhat_sq motheduc), yline(0) name(motheduc)
scatter (uhat_sq fatheduc), yline(0) name(enrol)
scatter (uhat_sq abil), yline(0) name(staff)
scatter (uhat_sq urban), yline(0) name(urban)

gr combine motheduc enrol staff urban

// Park test

quietly reg educ motheduc fatheduc abil urban
predict uhat, residual

gen uhat_sq = uhat^2
gen luhat_sq = log(uhat_sq)
gen lmotheduc = log(motheduc)
gen lfatheduc = log(fatheduc)
gen labil = log(abil)

reg luhat_sq lmotheduc lfatheduc labil urban


scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(4,.05)
scalar list  N R2 NR2 chi2


//Glesjar test

quietly reg educ motheduc fatheduc abil urban
predict uhat, residual

gen abs_uhat = abs(uhat)

reg abs_uhat motheduc fatheduc abil urban

scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(4,.05)
scalar list N R2 NR2 chi2

// Gold field Quandt test


*** choose sub-samples from the omotheducervations ***

scalar c = 1230/5
scalar sample_size = (1230 - c)/2
scalar list c sample_size

sort abil

reg educ motheduc fatheduc abil urban in 1/616
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

reg educ motheduc fatheduc abil urban in 615/1230
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)


scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)

scalar list F_stat pvalue Fcrit

// Breush-Pagan Godfrey test
 
quietly reg educ motheduc fatheduc abil urban
predict uhat, residual

gen uhat_sq = uhat^2
qui summarize uhat_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(1230-4)
scalar list RSS sigma_sq

gen p = uhat_sq/sigma_sq

reg p motheduc fatheduc abil urban

scalar ESS =  45.406184
scalar THETA = 0.5*ESS
scalar list THETA

scalar chi2=invchi2tail(4,.05)
scalar list THETA chi2


// White test

quietly regress educ motheduc fatheduc abil urban
predict uhat, residual

gen uhat_sq = uhat^2
gen motheduc_sq = motheduc^2
gen fatheduc_sq = fatheduc^2
gen abil_sq = abil^2


reg uhat_sq motheduc fatheduc abil urban motheduc_sq fatheduc_sq abil_sq 

scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(7,.05)
scalar list N R2 NR2 chi2

*** Autocorrelation ***

// graphical method

quietly regress educ motheduc fatheduc abil urban

gen time = _n
tsset time

predict uhat, residual
gen uhat_lag1 = L1.uhat
tsline uhat
twoway (scatter uhat_lag1 uhat)


// Runs test

quietly reg educ motheduc fatheduc abil urban
predict uhat, residual

scatter uhat time, connect(l) yline(0)

runtest uhat, threshold(0)


// Durbin Watson test

quietly regress educ motheduc fatheduc abil urban

estat dwatson

// Breusch-Godfrey test

quietly regress educ motheduc fatheduc abil urban

estat bgodfrey, lags(1/2)

*or

quietly regress educ motheduc fatheduc abil urban
predict uhat, residuals

gen uhat_lag1 = L1.uhat
gen uhat_lag2 = L2.uhat

reg uhat uhat_lag1 uhat_lag2, noconstant

scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2=invchi2tail(2,.05)
scalar list N R2 NR2 chi2


// Multicollinearity test

* looking at R-squared and t value

regress educ motheduc fatheduc abil urban


* Pair-wise correlations

pwcorr motheduc fatheduc abil urban

* Auxiliary regression for multicollinearity

// for motheduc
reg motheduc fatheduc abil urban

scalar R2 = _result(7)
scalar F_dist = (R2/(3-2))/((1-R2)/(3-2+1))
scalar F_crit = invFtail(3,3,.05)
scalar list R2 F_dist F_crit

// for fatheduc
reg fatheduc motheduc abil urban

scalar R2 = _result(7)
scalar F_dist = (R2/(3-2))/((1-R2)/(3-2+1))
scalar F_crit = invFtail(3,3,.05)
scalar list R2 F_dist F_crit

// for urban
reg urban  fatheduc motheduc abil 

scalar R2 = _result(7)
scalar F_dist = (R2/(3-2))/((1-R2)/(3-2+1))
scalar F_crit = invFtail(3,3,.05)
scalar list R2 F_dist F_crit

* partial correlations

pcorr motheduc fatheduc abil urban

* condition index

net install collin

collin motheduc fatheduc abil urban


* Ramsey reset test

quietly regress educ motheduc fatheduc abil urban

// Stata's version of the ramsey reset test

ovtest

// equivlent LM test


reg educ motheduc fatheduc abil urban

predict yhat, xb
gen yhat_2 = yhat^2
gen yhat_3 = yhat^3
gen yhat_4 = yhat^4

reg educ motheduc fatheduc abil urban yhat_2 yhat_3 yhat_4

scalar N=_result(1)
scalar R2_old= 0.4295
scalar R2_new = 0.4461 
scalar new_reg = 3
scalar no_param = 8
scalar NR2=N*R2
scalar F_value = ((R2_new - R2_old) / new_reg) / ((1 - R2_new) / (N - no_param))

scalar F_crit = invFtail(7,4,.05)

scalar P_value = Ftail(7,4,F_value)

scalar list R2_old R2_new N new_reg no_param F_value F_crit P_value