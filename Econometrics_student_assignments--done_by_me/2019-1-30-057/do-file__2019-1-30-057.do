use "link", clear

generate lbwght = log(bwght)
generate lfaminc = log(faminc)

reg lbwght cigs lfaminc parity male

** Assignment -- 1: graphical method

quietly reg lbwght cigs lfaminc parity male
predict uhat, resid

rvfplot, yline(0)

generate uhat_sq = uhat^2
scatter (uhat_sq cigs), yline(0) name(cigs)
scatter (uhat_sq lfaminc), yline(0) name(lfaminc)
scatter (uhat_sq parity), yline(0) name(parity)
scatter (uhat_sq male), yline(0) name(male)
graph combine cigs lfaminc parity male

** Assignment -- 2: Park test

quietly reg lbwght cigs lfaminc parity male
predict uhat, residual

generate uhat_sq = uhat^2
generate luhat_sq = log(uhat_sq)
generate lcigs = log(cigs)
generate lparity = log(parity)
generate llfaminc = log(lfaminc)

reg luhat_sq lcigs llfaminc lparity male


//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(4,.05)
scalar list N R2 LM crit_val

** Assignment -- 3:Glesjar test

quietly reg lbwght cigs lfaminc parity male
predict uhat, residual

generate abs_uhat = abs(uhat)
reg abs_uhat cigs lfaminc parity male

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(4,.05)
scalar list N R2 LM crit_val

** Assignment -- 4: Gold field Quandt test


// seclecting sub sample

scalar c = 1388/5
scalar sample_size = (1388 - c)/2
scalar list c sample_size

sort lfaminc

reg lbwght cigs lfaminc parity male in 1/556
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

reg lbwght cigs lfaminc parity male in 833/1388
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit

** Assignment -- 5: Breush-Pagan Godfrey test
 
quietly reg lbwght cigs lfaminc parity male
predict uhat, residual

generate uhat_sq = uhat^2
qui summarize uhat_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(1388-4)
scalar list RSS sigma_sq

generate p = uhat_sq/sigma_sq
reg p cigs lfaminc parity male

scalar ESS = 13.0231344
scalar THETA = 0.5*ESS

scalar crit_val=invchi2tail(4,.05)
scalar list THETA crit_val

// STATA - BG test
estat hettest

** Assignment -- 6: White test

quietly reg lbwght cigs lfaminc parity male
predict uhat, residual

generate uhat_sq = uhat^2
generate cigs_sq = cigs^2
generate lfaminc_sq = lfaminc^2
generate parity_sq = parity^2

reg uhat_sq cigs lfaminc parity male cigs_sq lfaminc_sq parity_sq

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(7,.05)
scalar list N R2 LM crit_val

// STATA - white test
estat imtest, white


*** Autocorrelation test ***

** Assignment -- 7: graphical method for autocorrelation

quietly reg lbwght cigs lfaminc parity male
predict uhat, residual

generate time = _n
tsset time

generate uhat_lag1 = L1.uhat
twoway (scatter uhat_lag1 uhat) (lfit uhat_lag1 uhat)


** Assignment -- 8: Runs test

quietly reg lbwght cigs lfaminc parity male
predict uhat, residual

scatter uhat time, connect(l) yline(0)

runtest uhat, threshold(0)


** Assignment -- 9: Durbin Watson test

quietly reg lbwght cigs lfaminc parity male
estat dwatson

** Assignment -- 10: Breusch-Godfrey test

quietly reg lbwght cigs lfaminc parity male
predict uhat, residuals

generate uhat_lag1 = L1.uhat
generate uhat_lag2 = L2.uhat

reg uhat uhat_lag1 uhat_lag2, noconstant

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2_crit=invchi2tail(2,.05)
scalar list N R2 LM chi2_crit

//In stata
quietly reg lbwght cigs lfaminc parity male
estat bgodfrey, lags(1/2)

*** Multicollinearity test ***

Assignment -- 11: looking at R-squared and t value

reg lbwght cigs lfaminc parity male


** Assignment -- 12: Pair-wise correlations

pwcorr cigs lfaminc parity male

** Assignment -- 13: Auxiliary reg for multicollinearity

// for cigs
reg cigs lfaminc parity male

// for lfaminc
reg lfaminc cigs parity male

// for male
reg male  lfaminc cigs parity 

// for parity
reg parity male  lfaminc cigs


** Assignment -- 14: partial correlations

pcorr cigs lfaminc parity male
pcorr lfaminc parity male cigs
pcorr parity cigs lfaminc male
pcorr male cigs lfaminc parity


** Assignment -- 15: condition index

collin cigs lfaminc parity male



** Assignment -- 16: Ramsey reset test

// Mannually
reg lbwght cigs lfaminc parity male

predict yhat, xb
generate yhat2 = yhat^2
generate yhat3 = yhat^3
generate yhat4=yhat^4

reg lbwght cigs lfaminc parity male yhat2 yhat3 yhat4

scalar N=_result(1)
scalar R2_old=  0.0350
scalar R2_new = 0.0396
scalar new_var = 3
scalar tot_param = 8

scalar F_value = ((R2_new - R2_old) / new_var) / ((1 - R2_new) / (N - tot_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)
scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value


// Stata's version of the ramsey reset test
quietly reg lbwght cigs lfaminc parity male
estat ovtest