clear all

import excel "C:\Users\nafew\Downloads\wti-remit-cpi (1).xlsx", sheet("Sheet1") firstrow

drop E F G H J K L M N O P Q S T U V W X Y Z I R

drop in 108/L

** set time variable
generate time = q(1993q3) + _n -1
format time %tq
tsset time

gen lcpi = log(cpi)
gen lremit = log(remit)
gen lwti = log(wti)

reg lcpi lremit lwti
estat sbsingle

tsline lcpi lremit lwti, name(A,replace)
twoway (line lcpi time, c(l) yaxis(1)) (line lremit time, c(l) yaxis(2)) (line lwti time, c(l) yaxis(3)), name(B,replace)
gr combine A B

** Step #1: DF-GLS unit root test
dfgls lcpi
dfgls lremit
dfgls lwti

** Step #2a: Engle-Granger Cointegration test
egranger lcpi lremit lwti, lags(1) 
// test for cointegration
egranger lcpi lremit lwti, lags(1) regress 
// LR estimates plus residual
egranger lcpi lremit lwti, lags(1) ecm 
// error-correction model

** Step #2: Johansen Cointegration test
varsoc lcpi lremit lwti
vecrank lcpi lremit lwti, lags(2) 
vec lcpi lremit lwti, lags(2)


** ans. 2 **

import excel "C:\Users\nafew\Downloads\DSE (1).xlsx", sheet("Sheet1") firstrow case(lower) clear
// to exclude weekends

bcal create calender, from(date) generate(b_date) maxgap(130) replace

tsset b_date


gen return = ln(dsex/l.dsex) * 100

gen return2 = return^2

tsline return
tsset newdate
format newdate %td

** ARCH(1)
arch return, arch(1)
est store ARCH
arch return, arch(1) vce(robust)
est store ARCH_rob
arch return, arch(1) vce(robust) distribution(t)
est store ARCH_t
arch return, arch(1) vce(robust) distribution(ged)
est store ARCH_ged
est table ARCH ARCH_rob ARCH_t ARCH_ged, se
est stats ARCH ARCH_rob ARCH_t ARCH_ged
** Student's t-distribution is best
arch return, arch(1) distribution(t) ar(1/5)
est store ARCH_arima
est stats ARCH_t ARCH_arima
** Student's t-distribution and ARIMA is best

** Day of the week effect
gen dayofweek = dow(newdate)
gen monday = 0
replace monday = 1 if dayofweek==1
gen friday = 0
replace friday = 1 if dayofweek==5
arch return, arch(1) distribution(t) ar(1/2) het(monday friday)
estat ic
** The model with Student's t-distribution and ARIMA is best

** Residual diagnostic checcking
*arch return, arch(1)
arch return, arch(1/10) // ARCH(p) model
predict residual, residuals
predict variance, variance
gen stdresidual = residual/variance^0.5
reg stdresidual
estat archlm, lags(1,5,10)
ac stdresidual
** There is ARCH effect. Simple ARCH(1) isn't adequate

** GACRH model 
arch return, arch(1) garch(1)
predict variance, variance
tsline variance
arch return, arch(1/5) garch(1/5)
arch return, arch(1) garch(1/2)






