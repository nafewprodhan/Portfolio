use "{{link}}", clear


regress price nox crime rooms stratio

*** Heteroskedasticity tesst ***

// A1: graphical method 

quietly regress price nox crime rooms stratio
predict err, resid

rvfplot, yline(0)

generate err_sq = err^2
scatter (err_sq nox), yline(0) name(nox)
scatter (err_sq crime), yline(0) name(crime)
scatter (err_sq rooms), yline(0) name(rooms)
scatter (err_sq stratio), yline(0) name(stratio)
graph combine nox crime rooms stratio

// A2: Park test

quietly regress price nox crime rooms stratio
predict err, residual

generate err_sq = err^2
generate lerr_sq = log(err_sq)
generate lnox = log(nox)
generate lrooms = log(rooms)
generate lcrime = log(crime)
generate lstratio = log(stratio)

regress lerr_sq lnox lcrime lrooms stratio

** LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi45=invchi2tail(4,.05)
scalar list  N R2 LM chi45


// A3: Glesjar test

quietly regress price nox crime rooms stratio
predict err, residual

generate abs_err = abs(err)

regress abs_err nox crime rooms stratio

** LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi45=invchi2tail(4,.05)
scalar list N R2 LM chi45

// A4: Gold field Quandt test

scalar c = 506/5
scalar sample_size = (506 - c)/2
scalar list c sample_size

sort nox

regress price nox crime rooms stratio in 1/203
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

regress price nox crime rooms stratio in 304/506
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)

scalar list F_stat pvalue Fcrit


// A5: Breush-Pagan Godfrey test
 
quietly regress price nox crime rooms stratio
predict err, residual

generate err_sq = err^2
qui summarize err_sq
scalar RSS = r(sum)
scalar sig_sq = RSS/(506-4)
scalar list RSS sig_sq

generate p_i = err_sq/sig_sq
regress p_i nox crime rooms stratio

scalar ESS = 108.723286 
scalar Theta = 0.5*ESS

scalar chi2=invchi2tail(4,.05)
scalar list Theta chi2

** BPG test in STATA
estat hettest

// A6: White test

quietly regress price nox crime rooms stratio
predict err, residual

generate err_sq = err^2
generate nox_sq = nox^2
generate crime_sq = crime^2
generate rooms_sq = rooms^2
generate stratio_sq = stratio^2

regress err_sq nox crime rooms stratio nox_sq crime_sq rooms_sq stratio_sq

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi75=invchi2tail(8,.05)
scalar list N R2 LM chi75

// Stata white test
estat imtest, white

*** Autocorrelation test ***

// A7: graphical method

quietly regress price nox crime rooms stratio
predict err, residual

generate time = _n
tsset time

generate err_lag1 = L1.err
twoway (scatter err_lag1 err) (lfit err_lag1 err)


// A8: Runs test

quietly regress price nox crime rooms stratio
predict err, residual

scatter err time, connect(l) yline(0)
runtest err, threshold(0)


// A9: Durbin Watson test

quietly regress price nox crime rooms stratio
estat dwatson

// A10: Breusch-Godfrey test

quietly regress price nox crime rooms stratio
predict err, residuals

generate err_lag1 = L1.err
generate err_lag2 = L2.err

regress err err_lag1 err_lag2, noconstant

** LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi25=invchi2tail(2,.05)
scalar list N R2 LM chi25

** Stata version
quietly regress price nox crime rooms stratio
estat bgodfrey, lags(1/2)

*** Multicollinearity test ***

// A11: looking at R-squared and t value

regress price nox crime rooms stratio


// A12: Pair-wise correlations

pwcorr nox crime rooms stratio

// A13: Auxiliary regress for multicollinearity

* nox
regress nox crime rooms stratio

* crime
regress crime nox rooms stratio

* stratio
regress stratio  crime nox rooms 

* rooms
regress rooms stratio  crime nox


// A14: partial correlations

pcorr nox crime rooms stratio
pcorr crime rooms stratio nox
pcorr rooms nox crime stratio
pcorr stratio nox crime rooms

// A15: condition index

net install collin
collin nox crime rooms stratio

// A16: Ramsey reset test




** Mannual test

regress price nox crime rooms stratio

predict yhat, xb
generate yhat2 = yhat^2
generate yhat3 = yhat^3
generate yhat4 = yhat^4

regress price nox crime rooms stratio yhat2 yhat3 yhat4

scalar N=_result(1)
scalar R2_old= 0.6140
scalar R2_new = 0.7041
scalar new_var = 3
scalar tot_param = 8
scalar F_value = ((R2_new - R2_old) / new_var) / ((1 - R2_new) / (N - tot_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)

scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value

// Stata's version
quietly regress price nox crime rooms stratio
estat ovtest


