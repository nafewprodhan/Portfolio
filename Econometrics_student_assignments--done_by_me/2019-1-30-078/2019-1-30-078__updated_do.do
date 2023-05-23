use "link", clear



generate lavgsal = log(avgsal)
generate lenrol = log(enrol)
generate lstaff = log(staff)

reg lavgsal bs lenrol lstaff lunch

// graphical method 

quietly reg lavgsal bs lenrol lstaff lunch
predict epsilon, resid

rvfplot, yline(0)

generate epsilon_sq = epsilon^2
scatter (epsilon_sq bs), yline(0) name(bs)
scatter (epsilon_sq lenrol), yline(0) name(lenrol)
scatter (epsilon_sq lstaff), yline(0) name(lstaff)
scatter (epsilon_sq lunch), yline(0) name(lunch)
graph combine bs lenrol lstaff lunch

// Park test
quietly reg lavgsal bs lenrol lstaff lunch
predict epsilon, residual

generate epsilon_sq = epsilon^2
generate lepsilon_sq = log(epsilon_sq)
generate lbs = log(bs)
generate llunch = log(lunch)
reg lepsilon_sq lbs lenrol lstaff lunch

** LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2=invchi2tail(4,.05)
scalar list  N R2 LM chi2


//Glesjar test

quietly reg lavgsal bs lenrol lstaff lunch
predict epsilon, residual

generate abs_epsilon = abs(epsilon)

reg abs_epsilon bs lenrol lstaff lunch

** LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2=invchi2tail(4,.05)
scalar list N R2 LM chi2

// Gold field Quandt test


*** choosing sub-samples ***

scalar c = 1848/5
scalar sample_size = (1848 - c)/2
scalar list c sample_size

sort bs

reg lavgsal bs lenrol lstaff lunch in 1/740
scalar rss1 = e(rmse)^2   
scalar df_rss1 = e(df_r)

reg lavgsal bs lenrol lstaff lunch in 1109/1848
scalar rss2 = e(rmse)^2
scalar df_rss2 = e(df_r)


scalar F_stat = rss2/rss1
scalar Fcrit = invFtail(df_rss2,df_rss1,.05)
scalar pvalue = Ftail(df_rss2,df_rss1,F_stat)
scalar list F_stat pvalue Fcrit

// Breush-Pagan Godfrey test
 
quietly reg lavgsal bs lenrol lstaff lunch
predict epsilon, residual

generate epsilon_sq = epsilon^2
qui summarize epsilon_sq
scalar RSS = r(sum)
scalar sigma_sq = RSS/(1848-4)
scalar list RSS sigma_sq

generate p = epsilon_sq/sigma_sq

reg p bs lenrol lstaff lunch

scalar ESS = 191.612285
scalar Theta = 0.5*ESS
scalar list Theta

scalar chi2=invchi2tail(4,.05)
scalar list Theta chi2

** STATA version BPG test for Heteroskedasticity
estat hettest

// White test

quietly regress lavgsal bs lenrol lstaff lunch
predict epsilon, residual

generate epsilon_sq = epsilon^2
generate bs_sq = bs^2
generate lenrol_sq = lenrol^2
generate lstaff_sq = lstaff^2
generate lunch_sq = lunch^2

reg epsilon_sq bs lenrol lstaff lunch bs_sq lenrol_sq lstaff_sq lunch_sq

** LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2=invchi2tail(8,.05)
scalar list N R2 LM chi2

** or
estat imtest, white


*** Autocorrelation ***

// graphical method to detect Autocorrelation

quietly regress lavgsal bs lenrol lstaff lunch

generate time = _n
tsset time

predict epsilon, residual
generate epsilon_lag1 = L1.epsilon
twoway (scatter epsilon_lag1 epsilon) (lfit epsilon_lag1 epsilon)


// Runs test

quietly reg lavgsal bs lenrol lstaff lunch
predict epsilon, residual

scatter epsilon time, connect(l) yline(0)
runtest epsilon, threshold(0)


// Durbin Watson test

quietly regress lavgsal bs lenrol lstaff lunch
estat dwatson

// Breusch-Godfrey test

quietly regress lavgsal bs lenrol lstaff lunch
estat bgodfrey, lags(1/2)

*or

quietly regress lavgsal bs lenrol lstaff lunch
predict epsilon, residuals

generate epsilon_lag1 = L1.epsilon
generate epsilon_lag2 = L2.epsilon

reg epsilon epsilon_lag1 epsilon_lag2, noconstant

** LM test
scalar N=_result(1)
scalar R2=_result(7)
scalar LM=N*R2
scalar chi2=invchi2tail(2,.05)
scalar list N R2 LM chi2


// Multicollinearity test

* By looking at R-squared and t value

reg lavgsal bs lenrol lstaff lunch


* Pair-wise correlations

pwcorr bs lenrol lstaff lunch

* Auxiliary regression for multicollinearity

// for bs
reg bs lenrol lstaff lunch

// for lenrol
reg lenrol bs lstaff lunch

// for lunch
reg lunch  lenrol bs lstaff 

// for lstaff
reg lstaff lenrol bs lunch

* partial correlations

pcorr bs lenrol lstaff lunch
pcorr lenrol lstaff lunch bs
pcorr lstaff bs lenrol lunch
pcorr lunch bs lenrol lstaff

* condition index

collin bs lenrol lstaff lunch


* Ramsey reset test

quietly regress lavgsal bs lenrol lstaff lunch

// Stata's version of the ramsey reset test

estat ovtest

// Mannual test

reg lavgsal bs lenrol lstaff lunch

predict yhat, xb
generate yhat2 = yhat^2
generate yhat3 = yhat^3
generate yhat4 = yhat^4

reg lavgsal bs lenrol lstaff lunch yhat2 yhat3 yhat4

scalar N=_result(1)
scalar R2_old= 0.4895
scalar R2_new =  0.4902 
scalar new_var = 3
scalar tot_param = 8
scalar F_value = ((R2_new - R2_old) / new_reg) / ((1 - R2_new) / (N - no_param))
scalar F_crit = invFtail(6,4,.05)
scalar P_value = Ftail(6,4,F_value)

scalar list R2_old R2_new N new_var tot_param F_value F_crit P_value