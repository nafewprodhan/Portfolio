use "link\2019-2-30-059.DTA", clear

regress cumgpa hsperc sat tothrs female

*** A1: graphical method 

quietly regress cumgpa hsperc sat tothrs female
predict error, resid

rvfplot, yline(0)

generate error_sq = error^2
scatter (error_sq hsperc), yline(0) name(hsperc)
scatter (error_sq sat), yline(0) name(sat)
scatter (error_sq tothrs), yline(0) name(tothrs)
scatter (error_sq female), yline(0) name(female)
graph combine hsperc sat tothrs female

*** A2: Park test

quietly regress cumgpa hsperc sat tothrs female
predict error, residual

generate error_sq = error^2
generate lerror_sq = log(error_sq)
generate lhsperc = log(hsperc)
generate ltothrs = log(tothrs)
generate lsat = log(sat)

regress lerror_sq lhsperc lsat ltothrs female

// LM for Park test 
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi45=invchi2tail(4,.05)
scalar list  N R2 LM chi2


*** A3: Glesjer test

quietly regress cumgpa hsperc sat tothrs female
predict error, residual

generate abs_error = abs(error)
regress abs_error hsperc sat tothrs female

// LM for Glesjer test 
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi45=invchi2tail(4,.05)
scalar list N R2 LM chi2


*** A4:  Gold field Quandt test


// choosing sub-samples

scalar c = 732/5
scalar sample_size = (732 - c)/2
scalar list c sample_size

sort sat

regress cumgpa hsperc sat tothrs female in 1/293
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

regress cumgpa hsperc sat tothrs female in 440/732
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit


*** A5:  Breush-Pagan Godfrey test
 
quietly regress cumgpa hsperc sat tothrs female
predict error, residual

generate error_sq = error^2
qui summarize error_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(732-4)
scalar list RSS sigma_sq

generate p = error_sq/sigma_sq

regress p hsperc sat tothrs female

scalar ESS = 598.233303 
scalar Theta = 0.5*ESS

scalar chi2=invchi2tail(4,.05)
scalar list Theta chi2

// stata version BP-godfrey test
quietly regress cumgpa hsperc sat tothrs female
estat hettest


*** A6: White test

quietly regress cumgpa hsperc sat tothrs female
predict error, residual

generate error_sq = error^2
generate hsperc_sq = hsperc^2
generate sat_sq = sat^2
generate tothrs_sq = tothrs^2

regress error_sq hsperc sat tothrs female hsperc_sq sat_sq tothrs_sq

//LM for white test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi75=invchi2tail(7,.05)
scalar list N R2 LM chi75

**** Autocorrelation ****

***A7: Graphical method

quietly regress cumgpa hsperc sat tothrs female
predict error, residual

generate time = _n
tsset time

generate error_lag1 = L1.error
twoway (scatter error_lag1 error) (lfit error_lag1 error)


***A8:  Runs test

quietly regress cumgpa hsperc sat tothrs female
predict error, residual

scatter error time, connect(l) yline(0) // to see the plot
runtest error, threshold(0)


***A9: Durbin Watson test -- By stata

quietly regress cumgpa hsperc sat tothrs female
estat dwatson

***A10:  Breusch-Godfrey test

//STATA version
quietly regress cumgpa hsperc sat tothrs female
estat bgodfrey, lags(1/2)

*or

//Mannual way of doing Breusch-Godfrey test
quietly regress cumgpa hsperc sat tothrs female
predict error, residuals

generate error_lag1 = L1.error
generate error_lag2 = L2.error

regress error error_lag1 error_lag2

// LM test for Breusch-Godfrey AC test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi25=invchi2tail(2,.05)
scalar list N R2 LM chi25


**** Multicollinearity test ****

***A11: looking at R-squared and t value

regress cumgpa hsperc sat tothrs female


***A12: Pair-wise correlations

pwcorr hsperc sat tothrs female


***A13: Auxiliary regress for multicollinearity

// for hsperc
regress hsperc sat tothrs female

// for sat
regress sat hsperc tothrs female

// for female
regress female  sat hsperc tothrs 

// for tothrs
regress tothrs female  sat hsperc  


***A10: partial correlations

pcorr hsperc sat tothrs female
pcorr sat tothrs female hsperc
pcorr tothrs hsperc sat female
pcorr female hsperc sat tothrs

***A10: condition index

collin hsperc sat tothrs female


* Ramsey reset test

// Stata's version of the ramsey reset test

quietly regress cumgpa hsperc sat tothrs female
ovtest


// equivlent mannual test
regress cumgpa hsperc sat tothrs female

predict yhat, xb
generate yhat2 = yhat^2
generate yhat3 = yhat^3
generate yhat4 = yhat^4

regress cumgpa hsperc sat tothrs female yhat2 yhat3 yhat4

scalar N=_result(1)
scalar R2_old= 0.2405
scalar R2_new = 0.2559
scalar new_var = 3
scalar tot_param = 8

scalar F_value = ((R2_new - R2_old) / new_reg) / ((1 - R2_new) / (N - tot_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)

scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value