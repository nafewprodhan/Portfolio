use "C:\Users\nafew\OneDrive - East West University\Academic\1. Fall-2022\A-465\Oboni\2019-3-30-029.DTA", clear

reg educ motheduc fatheduc abil wage

predict error, resid

* Graphical Method *

scatter (error educ), yline(0)

scatter (error motheduc), yline(0) name(motheduc)
scatter (error fatheduc), yline(0) name(fatheduc)
scatter (error abil), yline(0) name(abil)
scatter (error wage), yline(0) name(wage)

gr combine motheduc fatheduc abil wage


* Park test *

/*

h0: Homoskedasticity
ha: Heteroskedasticity


*/

reg educ motheduc fatheduc abil wage
predict error, resid

gen error_sq = error ^2

gen lerr_sq = log(error_sq)
gen lmotheduc = log(motheduc)
gen lfatheduc = log(fatheduc)
gen labil = log(abil)
gen lwage = log(wage)

reg lerr_sq lmotheduc lfatheduc labil lwage


* Glesjar test *

reg educ motheduc fatheduc abil wage
predict error, resid

gen abs_err = abs(error)

reg abs_err motheduc fatheduc abil wage


* Breusch Pagan test *

reg educ motheduc fatheduc abil wage
predict error, resid
gen error_sq = error ^2

reg error_sq motheduc fatheduc abil wage

scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar list N R2 NR2

scalar chi2=invchi2(4,.05)
scalar list chi2

// or

reg educ motheduc fatheduc abil wage

estat hettest


* white test *

reg educ motheduc fatheduc abil sibs
predict error, resid

gen error_sq = error ^2
gen medu_sq = motheduc^2
gen fedu_sq = fatheduc^2
gen abil_sq = abil^2
gen wage_sq = wage^2

reg error_sq motheduc medu_sq fatheduc fedu_sq abil abil_sq wage wage_sq

// extra

scalar N=_result(1)
scalar R2=_result(7)
scalar NR2=N*R2
scalar list N R2 NR2

scalar chi2=invchi2(8,.05)
scalar list chi2

//or

reg educ motheduc fatheduc abil wage

estat imtest, white


* Gold field quandt test *

sort abil

// 1/3 of total observation

reg educ motheduc fatheduc abil wage in 1/410
scalar rss1 = e(rmse)^2
scalar df_rss1 = e(df_r)

reg educ motheduc fatheduc abil wage in 820/1230
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

//F-statistics

scalar fgq = rss2/rss1

// F-critical

scalar fcrit = invFtail(df_rss2, df_rss1, .05)

// p-value

scalar pvalue = Ftail(df_rss2, df_rss1, fgq)

scalar list fgq pvalue fcrit

// If F-statistics > F-critical ==> presence of heteroskedasticity


**** AUTOCORRELATION ****

reg educ motheduc fatheduc abil sibs

predict error, resid

gen time = _n
tsset time

gen error_lag1 = L1.error

scatter error_lag1 error

// durbin watson test

reg educ motheduc fatheduc abil sibs

estat dwatson 


// B-godfrey test

reg educ motheduc fatheduc abil sibs

estat bgodfrey


// Runtest

reg educ motheduc fatheduc abil sibs
predict error, resid

scatter error time, connect(l) yline(0)

runtest error, threshold(0)












































































