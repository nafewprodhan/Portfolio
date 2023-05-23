use "link\2019-1-30-058.DTA", clear

regress bwght cigs fmaps npvis male

*** Assignment: 1 - graphical method

quietly regress bwght cigs fmaps npvis male
predict residual, resid

rvfplot, yline(0)

gen residual_sq = residual^2
scatter (residual_sq cigs), yline(0) name(cigs)
scatter (residual_sq fmaps), yline(0) name(fmaps)
scatter (residual_sq npvis), yline(0) name(npvis)
scatter (residual_sq male), yline(0) name(male)
graph combine cigs fmaps npvis male


*** Assignment: 2 - Park test

quietly regress bwght cigs fmaps npvis male
predict residual, residual

gen residual_sq = residual^2
gen lresidual_sq = log(residual_sq)
gen lcigs = log(cigs)
gen lnpvis = log(npvis)
gen lfmaps = log(fmaps)

regress lresidual_sq lcigs lfmaps lnpvis male

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar chi2_crit=invchi2tail(4,.05)
scalar list N R2 NR2 chi2_crit

*** Assignment: 3 - Glesjar test

quietly regress bwght cigs fmaps npvis male
predict residual, residual

gen abs_residual = abs(residual)
regress abs_residual cigs fmaps npvis male

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar CRIT_VAL=invchi2tail(4,.05)
scalar list N R2 NR2 CRIT_VAL

*** Assignment: 4 -  Gold field Quandt test


// to select the sub samles

drop if residual == .

scalar c = 1653/5
scalar sample_size = (1653 - c)/2
scalar list c sample_size

// sorting
sort fmaps

// generating RSS/df
regress bwght cigs fmaps npvis male in 1/662
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

regress bwght cigs fmaps npvis male in 992/1653
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

// comaring F-stat to F-crit
scalar F_stat = rss2/rss1
scalar F_crit = invFtail(df_rss2,df_rss1,.05)
scalar P_value = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat P_value F_crit

*** Assignment: 5 -  Breush-Pagan Godfrey test
 
quietly regress bwght cigs fmaps npvis male
predict residual, residual

gen residual_sq = residual^2
qui summarize residual_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(1832-4)
scalar list RSS sigma_sq

gen p = residual_sq/sigma_sq
regress p cigs fmaps npvis male

scalar ESS = 315.104322
scalar THETA = 0.5*ESS

scalar CRIT_VAL=invchi2tail(4,.05)
scalar list THETA CRIT_VAL

// stata BPG test
estat hettest

*** Assignment: 6 -  White test

quietly regress bwght cigs fmaps npvis male
predict residual, residual

gen residual_sq = residual^2
gen cigs_sq = cigs^2
gen fmaps_sq = fmaps^2
gen npvis_sq = npvis^2

regress residual_sq cigs fmaps npvis male cigs_sq fmaps_sq npvis_sq

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar CRIT_VAL=invchi2tail(7,.05)
scalar list N R2 NR2 CRIT_VAL

// Stata white test
estat imtest, white

*****  Autocorrelation Test  *****

*** Assignment: 7 -  graphical method for AC

quietly regress bwght cigs fmaps npvis male
predict residual, residual

gen time = _n
tsset time

gen residual_lag1 = L1.residual
twoway (scatter residual_lag1 residual) (lfit residual_lag1 residual)


*** Assignment: 8 -  Runs test

quietly regress bwght cigs fmaps npvis male
predict residual, residual

scatter residual time, connect(l) yline(0)
runtest residual, threshold(0)


*** Assignment: 9 -  Durbin Watson test by stata

quietly regress bwght cigs fmaps npvis male
estat dwatson

*** Assignment: 10 -  Breusch-Godfrey test

//by stata
quietly regress bwght cigs fmaps npvis male
estat bgodfrey, lags(1/2)

//Mannual
quietly regress bwght cigs fmaps npvis male
predict residual, residuals

gen residual_lag1 = L1.residual
gen residual_lag2 = L2.residual

regress residual residual_lag1 residual_lag2, noconstant

scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar CRIT_VAL=invchi2tail(2,.05)
scalar list N R2 NR2 CRIT_VAL


***** Multicollinearity test ******

*** Assignment: 11 -  looking at R-squared and t value

regress bwght cigs fmaps npvis male


*** Assignment: 12 -  Pair-wise correlations
pwcorr cigs fmaps npvis male

*** Assignment: 13 -  Auxiliary regress for multicollinearity
// for cigs
regress cigs fmaps npvis male
// for fmaps
regress fmaps cigs npvis male
// for male
regress male  fmaps cigs npvis 
// for npvis
regress npvis male  fmaps cigs

*** Assignment: 14 -  partial correlations

pcorr cigs fmaps npvis male
pcorr fmaps npvis male cigs
pcorr npvis cigs fmaps male
pcorr male cigs fmaps npvis

*** Assignment: 15 -  condition index

net install collin
collin cigs fmaps npvis male


*** Assignment: 16 -  Ramsey reset test

//by stata
quietly regress bwght cigs fmaps npvis male
estat ovtest

// equivlent Mannual test
regress bwght cigs fmaps npvis male

predict yhat, xb
gen yhat2 = yhat^2
gen yhat3 = yhat^3
gen yhat4=yhat^4

regress bwght cigs fmaps npvis male yhat2 yhat3 yhat4

scalar N=_result(1)
scalar R2_old= 0.0466
scalar R2_new = 0.0605
scalar new_var = 3
scalar tot_param = 8

scalar F_value = ((R2_new - R2_old) / new_var) / ((1 - R2_new) / (N - tot_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)
scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value