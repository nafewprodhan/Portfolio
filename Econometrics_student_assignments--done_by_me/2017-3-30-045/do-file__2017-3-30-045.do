use "link", clear

reg crimes unem crimes_lag1 lawexpc nrtheast

**** Heteroskedasticity Test ****

** A:1 -- graphical method 

quietly reg crimes unem crimes_lag1 lawexpc nrtheast
predict err, resid

rvfplot, yline(0)

gen err_sq = err^2
scatter (err_sq unem), yline(0) name(unem)
scatter (err_sq crimes_lag1), yline(0) name(crimes_lag1)
scatter (err_sq lawexpc), yline(0) name(lawexpc)
scatter (err_sq nrtheast), yline(0) name(nrtheast)
graph combine unem crimes_lag1 lawexpc nrtheast

** A:2 -- Park test

quietly reg crimes unem crimes_lag1 lawexpc nrtheast
predict err, residual

gen err_sq = err^2
gen lerr_sq = log(err_sq)
gen lunem = log(unem)
gen llawexpc = log(lawexpc)
gen lcrimes_lag1 = log(crimes_lag1)

reg lerr_sq lunem lcrimes_lag1 llawexpc nrtheast

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(4,.05)
scalar list N R2 LM crit_val

** A:3 -- Glesjar test

quietly reg crimes unem crimes_lag1 lawexpc nrtheast
predict err, residual

gen abs_err = abs(err)
reg abs_err unem crimes_lag1 lawexpc nrtheast

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(4,.05)
scalar list N R2 LM crit_val

** A:4 --  Gold field Quandt test


// to choose sample size ***

scalar c = 92/5
scalar sample_size = (92 - c)/2
scalar list c sample_size

sort unem

reg crimes unem crimes_lag1 lawexpc nrtheast in 1/37
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

reg crimes unem crimes_lag1 lawexpc nrtheast in 55/92
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)

scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit

** A:5 -- Breush-Pagan Godfrey test
 
quietly reg crimes unem crimes_lag1 lawexpc nrtheast
predict err, residual

gen err_sq = err^2
qui summarize err_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(92-4)
scalar list RSS sigma_sq

gen p = err_sq/sigma_sq
reg p unem crimes_lag1 lawexpc nrtheast

scalar ess = 76.9850646
scalar theta = 0.5*ESS

scalar crit_val=invchi2tail(4,.05)
scalar list theta crit_val

// Stata version
estat hettest


** A:6 --  White test

quietly reg crimes unem crimes_lag1 lawexpc nrtheast
predict err, residual

gen err_sq = err^2
gen unem_sq = unem^2
gen crimes_lag1_sq = crimes_lag1^2
gen lawexpc_sq = lawexpc^2

reg err_sq unem crimes_lag1 lawexpc nrtheast unem_sq crimes_lag1_sq lawexpc_sq

// LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar crit_val=invchi2tail(7,.05)
scalar list N R2 LM crit_val

// stata version of white test
estat imtest, white


*** Autocorrelation test ***

** A:7 --  graphical method

quietly reg crimes unem crimes_lag1 lawexpc nrtheast
predict err, residual

gen time = _n
tsset time


gen err_lag1 = L1.err
twoway (scatter err_lag1 err) (lfit err_lag1 err)


** A:8 --  Runs test

quietly reg crimes unem crimes_lag1 lawexpc nrtheast
predict err, residual

scatter err time, connect(l) yline(0)
runtest err, threshold(0)


** A:9 --  Durbin Watson test

//STATA version
quietly reg crimes unem crimes_lag1 lawexpc nrtheast
estat dwatson

** A:10 --  Breusch-Godfrey test

//STATA version
quietly reg crimes unem crimes_lag1 lawexpc nrtheast
estat bgodfrey, lags(1/2)

*or mannual BG test

quietly reg crimes unem crimes_lag1 lawexpc nrtheast
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


***** Multicollinearity test *****

** A:11 --  By looking at R-squared and t value

reg crimes unem crimes_lag1 lawexpc nrtheast


** A:12 --  Pair-wise correlations

pwcorr unem crimes_lag1 lawexpc nrtheast

** A:13 --  Auxiliary reg for multicollinearity

// for unem
reg unem crimes_lag1 lawexpc nrtheast

// for crimes_lag1
reg crimes_lag1 unem lawexpc nrtheast

// for nrtheast
reg nrtheast  crimes_lag1 unem lawexpc 

// for lawexpc
reg lawexpc nrtheast  crimes_lag1 unem


** A:14 --  partial correlations

pcorr unem crimes_lag1 lawexpc nrtheast
pcorr crimes_lag1 lawexpc nrtheast unem
pcorr lawexpc unem crimes_lag1 nrtheast
pcorr nrtheast unem crimes_lag1 lawexpc

** A:15 --  condition index

net install collin
collin unem crimes_lag1 lawexpc nrtheast


** A:16 --  Ramsey reset test

// Stata's version of the ramsey reset test
quietly reg crimes unem crimes_lag1 lawexpc nrtheast
estat ovtest

// equivlent Mannual Reset test

reg crimes unem crimes_lag1 lawexpc nrtheast

predict yhat, xb
gen yhat2 = yhat^2
gen yhat3 = yhat^3
gen yhat4=yhat^4

reg crimes unem crimes_lag1 lawexpc nrtheast yhat2 yhat3 yhat4

scalar N=_result(1)
scalar R2_old= 0.3716
scalar R2_new =  0.3802
scalar new_var = 3
scalar tot_param = 8

scalar F_value = ((R2_new - R2_old) / new_var) / ((1 - R2_new) / (N - tot_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)

scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value