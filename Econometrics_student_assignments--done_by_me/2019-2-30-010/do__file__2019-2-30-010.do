use "C:\Users\USER\OneDrive - East West University\Academic\1. Fall-2022\A-465\0. done\2019-2-30-010\2019-2-30-010.DTA", clear

reg cumgpa hsperc sat tothrs female

****** Heteroskedasticity test *******

** 1. graphical method 

quietly reg cumgpa hsperc sat tothrs female
predict err, resid

rvfplot, yline(0)

gen err_sq = err^2
scatter (err_sq hsperc), yline(0) name(hsperc)
scatter (err_sq sat), yline(0) name(sat)
scatter (err_sq tothrs), yline(0) name(tothrs)
scatter (err_sq female), yline(0) name(female)
graph combine hsperc sat tothrs female

** 2. Park test

quietly reg cumgpa hsperc sat tothrs female
predict err, residual

gen err_sq = err^2
gen lerr_sq = log(err_sq)
gen lhsperc = log(hsperc)
gen ltothrs = log(tothrs)
gen lsat = log(sat)

reg lerr_sq lhsperc lsat ltothrs female


//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2_crit=invchi2tail(4,.05)
scalar list  N R2 LM chi2_crit


** 3. Glesjer test

quietly reg cumgpa hsperc sat tothrs female
predict err, residual

gen abs_err = abs(err)

reg abs_err hsperc sat tothrs female

scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2_crit=invchi2tail(4,.05)
scalar list N R2 LM chi2_crit

** 4. Gold field Quandt test


// finding sub-samples

scalar c = 732/5
scalar sample_size = (732 - c)/2
scalar list c sample_size

sort sat

reg cumgpa hsperc sat tothrs female in 1/293
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

reg cumgpa hsperc sat tothrs female in 440/732
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)


scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)

scalar list F_stat pvalue Fcrit

** 5. Breush-Pagan Godfrey test
 
quietly reg cumgpa hsperc sat tothrs female
predict err, residual

gen err_sq = err^2
qui summarize err_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(732-4)
scalar list RSS sigma_sq

gen pi = err_sq/sigma_sq

reg pi hsperc sat tothrs female

scalar ESS = 598.233303 
scalar Theta = 0.5*ESS

scalar chi2=invchi2tail(4,.05)
scalar list Theta chi2

// Stata version of BG test
estat hettest

** 6. White test

quietly reg cumgpa hsperc sat tothrs female
predict err, residual

gen err_sq = err^2
gen hsperc_sq = hsperc^2
gen sat_sq = sat^2
gen tothrs_sq = tothrs^2

reg err_sq hsperc sat tothrs female hsperc_sq sat_sq tothrs_sq


// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2_crit=invchi2tail(7,.05)
scalar list N R2 LM chi2_crit

****** Autocorrelation test ******

** 7. Graphical method

quietly reg cumgpa hsperc sat tothrs female
predict err, residual

gen time = _n
tsset time

gen err_lag1 = L1.err
twoway (scatter err_lag1 err) (lfit err_lag1 err)


** 8. Runs test

quietly reg cumgpa hsperc sat tothrs female
predict err, residual

scatter err time, connect(l) yline(0)
runtest err, threshold(0)


** 9. Durbin Watson test

quietly reg cumgpa hsperc sat tothrs female
estat dwatson

** 10. Breusch-Godfrey test

//STATA version of BG test for auto correlation
quietly reg cumgpa hsperc sat tothrs female
estat bgodfrey, lags(1/2)

// Mannually

quietly reg cumgpa hsperc sat tothrs female
predict err, residuals

gen err_lag1 = L1.err
gen err_lag2 = L2.err

reg err err_lag1 err_lag2, noconstant

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2_crit=invchi2tail(2,.05)
scalar list N R2 LM chi2_crit


****** Multicollinearity test *******

** 11.By observing the value of R-square and t

reg cumgpa hsperc sat tothrs female


** 12. Pair-wise correlations

pwcorr hsperc sat tothrs female

** 13. Auxiliary reg for multicollinearity

reg hsperc sat tothrs female
reg sat hsperc tothrs female
reg female  sat hsperc tothrs
reg tothrs female  sat hsperc

** 14. partial correlations

pcorr hsperc sat tothrs female
pcorr sat tothrs female hsperc
pcorr tothrs hsperc sat female
pcorr female hsperc sat tothrs

** 15. condition index

net install collin
collin hsperc sat tothrs female


** 16.  Ramsey reset test

// Stata's version of the ramsey reset test

quietly reg cumgpa hsperc sat tothrs female
estat ovtest

// equivlent Manual test

reg cumgpa hsperc sat tothrs female

predict yhat, xb
gen yhat2 = yhat^2
gen yhat3 = yhat^3

reg cumgpa hsperc sat tothrs female yhat2 yhat3

scalar N=_result(1)
scalar R2_old= 0.2405
scalar R2_new = 0.2471
scalar new_var = 2
scalar tot_param = 7

scalar F_value = ((R2_new - R2_old) / new_var) / ((1 - R2_new) / (N - tot_param))
scalar F_crit = invFtail(6,4,.05)
scalar P_value = Ftail(6,4,F_value)

scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value