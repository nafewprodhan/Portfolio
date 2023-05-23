use "C:\Users\nafew\OneDrive - East West University\Academic\1. Fall-2022\A-465\2019-2-30-025\2019-2-30-025.DTA", clear


regress educ motheduc fatheduc abil urban

*** Assignment 1: Graphical Method for hettest

quietly regress educ motheduc fatheduc abil urban
predict epsilon, resid

rvfplot, yline(0)

generate epsilon_sq = epsilon^2
scatter (epsilon_sq motheduc), yline(0) name(motheduc)
scatter (epsilon_sq fatheduc), yline(0) name(enrol)
scatter (epsilon_sq abil), yline(0) name(staff)
scatter (epsilon_sq urban), yline(0) name(urban)
gr combine motheduc enrol staff urban

*** Assignment 2:  Park test

quietly regress educ motheduc fatheduc abil urban
predict epsilon, residual

generate epsilon_sq = epsilon^2
generate lepsilon_sq = log(epsilon_sq)
generate lmotheduc = log(motheduc)
generate lfatheduc = log(fatheduc)
generate labil = log(abil)
regress lepsilon_sq lmotheduc lfatheduc labil urban

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM_stat=N*R2
scalar CHI2_crit=invchi2tail(4,.05)
scalar list  N R2 LM_stat CHI2_crit


*** Assignment 3: Glesjar test

quietly regress educ motheduc fatheduc abil urban
predict epsilon, residual
generate abs_epsilon = abs(epsilon)

regress abs_epsilon motheduc fatheduc abil urban

//LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM_stat=N*R2
scalar CHI2_crit=invchi2tail(4,.05)
scalar list  N R2 LM_stat CHI2_crit

*** Assignment 4:  Gold field Quandt test

scalar c = 1230/5
scalar sample_size = (1230 - c)/2
scalar list c sample_size

sort abil

regress educ motheduc fatheduc abil urban in 1/616
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

regress educ motheduc fatheduc abil urban in 615/1230
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)


scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit

*** Assignment 5:  Breush-Pagan Godfrey test
 
quietly regress educ motheduc fatheduc abil urban
predict epsilon, residual

generate epsilon_sq = epsilon^2
quietly summarize epsilon_sq
scalar RSS = r(sum)
scalar SIGMA_SQ = RSS/(1230-4)
scalar list RSS SIGMA_SQ

generate pi = epsilon_sq/SIGMA_SQ
regress pi motheduc fatheduc abil urban

scalar ESS =  45.406184
scalar THETA = 0.5*ESS

scalar CHI2_crit=invchi2tail(4,.05)
scalar list THETA CHI2_crit

// STATA test for BPG
quietly regress educ motheduc fatheduc abil urban
estat hettest

*** Assignment 6:  White test general

quietly regress educ motheduc fatheduc abil urban
predict epsilon, residual

generate epsilon_sq = epsilon^2
generate motheduc_sq = motheduc^2
generate fatheduc_sq = fatheduc^2
generate abil_sq = abil^2

regress epsilon_sq motheduc fatheduc abil urban motheduc_sq fatheduc_sq abil_sq 

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

quietly regress educ motheduc fatheduc abil urban

generate time = _n
tsset time

predict epsilon, residual
generate epsilon_lag1 = L1.epsilon
tsline epsilon
twoway (scatter epsilon_lag1 epsilon)

*** Assignment 8:  Runs test

quietly regress educ motheduc fatheduc abil urban
predict epsilon, residual

scatter epsilon time, connect(l) yline(0)
runtest epsilon, threshold(0)


*** Assignment 9: Durbin Watson test

quietly regress educ motheduc fatheduc abil urban
estat dwatson

*** Assignment 10:  Breusch-Godfrey test

quietly regress educ motheduc fatheduc abil urban
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
quietly regress educ motheduc fatheduc abil urban
estat bgodfrey, lags(1/2)

*** Multicollinearity test ***

*** Assignment 11: looking at R-squared and t value

regress educ motheduc fatheduc abil urban


*** Assignment 12: Pair-wise correlations

pwcorr motheduc fatheduc abil urban

*** Assignment 13: Auxiliary regression for multicollinearity

// for motheduc
regress motheduc fatheduc abil urban
// for fatheduc
regress fatheduc motheduc abil urban
// for abil
regress abil urban  fatheduc motheduc
// for urban
regress urban  fatheduc motheduc abil 


*** Assignment 14: partial correlations

pcorr motheduc fatheduc abil urban
pcorr fatheduc motheduc  abil urban
pcorr urban fatheduc motheduc  abil 
pcorr abil fatheduc motheduc   urban

*** Assignment 15: condition index
net install collin
collin motheduc fatheduc abil urban


*** Assignment 16: Ramsey reset test


// equivlent mannual N-order Ramsey reset test

regress educ motheduc fatheduc abil urban

predict educ_hat, xb
generate educ_hat2 = educ_hat^2
generate educ_hat3 = educ_hat^3
generate educ_hat4 = educ_hat^4

regress educ motheduc fatheduc abil urban educ_hat2 educ_hat3 educ_hat4

scalar N=_result(1)
scalar R2_old= 0.4295
scalar R2_new = 0.4507 
scalar new_reg = 3
scalar no_param = 8
scalar NR2=N*R2

scalar F_value = ((R2_new - R2_old) / new_reg) / ((1 - R2_new) / (N - no_param))
scalar F_crit = invFtail(7,4,.05)
scalar P_value = Ftail(7,4,F_value)
scalar list R2_old R2_new N new_reg no_param F_value F_crit P_value

// Stata's version of the ramsey reset test
quietly regress educ motheduc fatheduc abil urban
estat ovtest