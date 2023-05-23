use "(link)\2019-1-30-015.DTA", clear

reg crmrte prbarr prbconv polpc west

*** A- 1: Graphical Method for hettest

quietly reg crmrte prbarr prbconv polpc west
predict error, resid

rvfplot, yline(0)


*** A- 2:  Park test

quietly reg crmrte prbarr prbconv polpc west
predict error, residual

gen error_sq = error^2
gen lerror_sq = log(error_sq)
gen lprbarr = log(prbarr)
gen lprbconv = log(prbconv)
gen lpolpc = log(polpc)
reg lerror_sq lprbarr lprbconv lpolpc west


*** A- 3: Glesjar test

quietly reg crmrte prbarr prbconv polpc west
predict error, residual
gen abs_error = abs(error)

reg abs_error prbarr prbconv polpc west


*** A- 4:  Gold field Quandt test

scalar c = 630/5
scalar sample_size = (630 - c)/2
scalar list c sample_size

sort prbconv

reg crmrte prbarr prbconv polpc west in 1/252
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

reg crmrte prbarr prbconv polpc west in 379/630
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)


scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit

*** A- 5:  Breush-Pagan Godfrey test
 
quietly reg crmrte prbarr prbconv polpc west
predict error, residual

gen error_sq = error^2
quietly summarize error_sq
scalar RSS = r(sum)
scalar SIGMA_SQ = RSS/(630-4)
scalar list RSS SIGMA_SQ

gen pi = error_sq/SIGMA_SQ
reg pi prbarr prbconv polpc west

scalar ESS =   478.898945 
scalar THETA = 0.5*ESS

scalar CHI2_crit=invchi2tail(4,.05)
scalar list THETA CHI2_crit

// STATA test for BPG
quietly reg crmrte prbarr prbconv polpc west
estat hettest

*** A- 6:  White test general

quietly reg crmrte prbarr prbconv polpc west
predict error, residual

gen error_sq = error^2
gen prbarr_sq = prbarr^2
gen prbconv_sq = prbconv^2
gen polpc_sq = polpc^2

reg error_sq prbarr prbconv polpc west prbarr_sq prbconv_sq polpc_sq 

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM_stat=N*R2
scalar CHI2_crit=invchi2tail(7,.05)
scalar list N R2 LM_stat CHI2_crit

// STATA white test
estat imtest, white


*** Autocorrelation Test ***

*** A- 7:  graphical method

quietly reg crmrte prbarr prbconv polpc west

gen time = _n
tsset time

predict error, residual
gen error_lag1 = L1.error
twoway (scatter error_lag1 error)

*** A- 8:  Runs test

quietly reg crmrte prbarr prbconv polpc west
predict error, residual

runtest error, threshold(0)


*** A- 9: Durbin Watson test

quietly reg crmrte prbarr prbconv polpc west
estat dwatson

*** A- 10:  Breusch-Godfrey test

quietly reg crmrte prbarr prbconv polpc west
predict error, residuals

gen error_lag1 = L1.error
gen error_lag2 = L2.error

reg error error_lag1 error_lag2, noconstant

scalar N=_result(1)
scalar R2=_result(7)
scalar LM_stat = N*R2
scalar CHI2_crit=invchi2tail(2,.05)
scalar list N R2 LM_stat CHI2_crit

//stata built in command for bg test
quietly reg crmrte prbarr prbconv polpc west
estat bgodfrey, lags(1)

*** Multicollinearity test ***

*** A- 11: looking at R-squared and t value

reg crmrte prbarr prbconv polpc west


*** A- 12: Pair-wise correlations

pwcorr prbarr prbconv polpc west

*** A- 13: Auxiliary regression for multicollinearity

// for prbarr
reg prbarr prbconv polpc west
// for prbconv
reg prbconv prbarr polpc west
// for polpc
reg polpc west  prbconv prbarr
// for west
reg west  prbconv prbarr polpc 


*** A- 14: partial correlations

pcorr prbarr prbconv polpc west
pcorr prbconv prbarr  polpc west
pcorr west prbconv prbarr  polpc 
pcorr polpc prbconv prbarr   west

*** A- 15: condition index

collin prbarr prbconv polpc west


*** A- 16: Ramsey reset test

reg crmrte prbarr prbconv polpc west

predict crmrte_hat, xb
gen crmrte_hat2 = crmrte_hat^2
gen crmrte_hat3 = crmrte_hat^3
gen crmrte_hat4 = crmrte_hat^4

reg crmrte prbarr prbconv polpc west crmrte_hat2 crmrte_hat3 crmrte_hat4

scalar N=_result(1)
scalar R2_old=  0.4026
scalar R2_new =  0.5391
scalar new_reg = 3
scalar no_param = 8

scalar F_value = ((R2_new - R2_old) / new_reg) / ((1 - R2_new) / (N - no_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)
scalar list R2_old R2_new N new_reg no_param F_value F_crit P_value

// Stata's version of the ramsey reset test
quietly reg crmrte prbarr prbconv polpc west
estat ovtest