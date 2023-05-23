use "C:\Users\nafew\OneDrive - East West University\Academic\1. Fall-2022\A-465\Tayan\2019-2-30-009.dta", clear

gen lwage = log(hrwage)

regress lwage educ earns exper hispanic

*** Assignment 1: Graphical Method for hettest

quietly regress lwage educ earns exper hispanic
predict epsilon, resid

rvfplot, yline(0)


*** Assignment 2:  Park test

quietly regress lwage educ earns exper hispanic
predict epsilon, residual

generate epsilon_sq = epsilon^2
generate lepsilon_sq = log(epsilon_sq)
generate leduc = log(educ)
generate learns = log(earns)
generate lexper = log(exper)
regress lepsilon_sq leduc learns lexper hispanic


*** Assignment 3: Glesjar test

quietly regress lwage educ earns exper hispanic
predict epsilon, residual
generate abs_epsilon = abs(epsilon)

regress abs_epsilon educ earns exper hispanic


*** Assignment 4:  Gold field Quandt test

scalar c = 3286/5
scalar sample_size = (3286 - c)/2
scalar list c sample_size

sort earns

regress lwage educ earns exper hispanic in 1/1315
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

regress lwage educ earns exper hispanic in 1972/3286
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)


scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit

*** Assignment 5:  Breush-Pagan Godfrey test
 
quietly regress lwage educ earns exper hispanic
predict epsilon, residual

generate epsilon_sq = epsilon^2
quietly summarize epsilon_sq
scalar RSS = r(sum)
scalar SIGMA_SQ = RSS/(3286-4)
scalar list RSS SIGMA_SQ

generate pi = epsilon_sq/SIGMA_SQ
regress pi educ earns exper hispanic

scalar ESS =   149.426885 
scalar THETA = 0.5*ESS

scalar CHI2_crit=invchi2tail(4,.05)
scalar list THETA CHI2_crit

// STATA test for BPG
quietly regress lwage educ earns exper hispanic
estat hettest

*** Assignment 6:  White test general

quietly regress lwage educ earns exper hispanic
predict epsilon, residual

generate epsilon_sq = epsilon^2
generate educ_sq = educ^2
generate earns_sq = earns^2
generate exper_sq = exper^2

regress epsilon_sq educ earns exper hispanic educ_sq earns_sq exper_sq 

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM_stat=N*R2
scalar CHI2_crit=invchi2tail(7,.05)
scalar list N R2 LM_stat CHI2_crit

// STATA white test
estat imtest, white


*** Autocorrelation Test ***

*** Assignment 7:  graphical method

quietly regress lwage educ earns exper hispanic

generate time = _n
tsset time

predict epsilon, residual
generate epsilon_lag1 = L1.epsilon
twoway (scatter epsilon_lag1 epsilon)

*** Assignment 8:  Runs test

quietly regress lwage educ earns exper hispanic
predict epsilon, residual

runtest epsilon, threshold(0)


*** Assignment 9: Durbin Watson test

quietly regress lwage educ earns exper hispanic
estat dwatson

*** Assignment 10:  Breusch-Godfrey test

quietly regress lwage educ earns exper hispanic
predict epsilon, residuals

generate epsilon_lag1 = L1.epsilon
generate epsilon_lag2 = L2.epsilon

regress epsilon epsilon_lag1 epsilon_lag2, noconstant

scalar N=_result(1)
scalar R2=_result(7)
scalar LM_stat = N*R2
scalar CHI2_crit=invchi2tail(2,.05)
scalar list N R2 LM_stat CHI2_crit

//stata built in command for bg test
quietly regress lwage educ earns exper hispanic
estat bgodfrey, lags(1/2)

*** Multicollinearity test ***

*** Assignment 11: looking at R-squared and t value

regress lwage educ earns exper hispanic


*** Assignment 12: Pair-wise correlations

pwcorr educ earns exper hispanic

*** Assignment 13: Auxiliary regression for multicollinearity

// for educ
regress educ earns exper hispanic
// for earns
regress earns educ exper hispanic
// for exper
regress exper hispanic  earns educ
// for hispanic
regress hispanic  earns educ exper 


*** Assignment 14: partial correlations

pcorr educ earns exper hispanic
pcorr earns educ  exper hispanic
pcorr hispanic earns educ  exper 
pcorr exper earns educ   hispanic

*** Assignment 15: condition index
net install collin
collin educ earns exper hispanic


*** Assignment 16: Ramsey reset test

regress lwage educ earns exper hispanic

predict lwage_hat, xb
generate lwage_hat2 = lwage_hat^2
generate lwage_hat3 = lwage_hat^3
generate lwage_hat4 = lwage_hat^4

regress lwage educ earns exper hispanic lwage_hat2 lwage_hat3 lwage_hat4

scalar N=_result(1)
scalar R2_old= 0.6593
scalar R2_new = 0.7017
scalar new_reg = 3
scalar no_param = 8

scalar F_value = ((R2_new - R2_old) / new_reg) / ((1 - R2_new) / (N - no_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)
scalar list R2_old R2_new N new_reg no_param F_value F_crit P_value

// Stata's version of the ramsey reset test
quietly regress lwage educ earns exper hispanic
estat ovtest