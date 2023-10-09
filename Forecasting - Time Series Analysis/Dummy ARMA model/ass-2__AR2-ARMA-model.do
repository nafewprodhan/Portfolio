// Modify the Stata code below to estimate an AR(2) and ARMA(2,1) model.

*******************************************
** Simulate AR(1) process
// clear all
// quietly set obs 100
// set seed 10101
// gen t=_n
// scalar a = 0.48
// scalar b = 0.83
// scalar c = 0.005
// scalar d = 0.014
// gen rho = a+(b-a)*runiform()
// gen s = c+(d-c)*runiform()
// gen et = rnormal(0, s)
// quietly gen yt= et in 1
// quietly replace yt=rho*yt[_n-1]+et in 2/L // rho = psi_1
// tsset t
// ac yt // indicates MA(3)
// pac yt // indicates AR(1)
* so an ARMA(1,3) model can be estimated
*******************************************

clear all
quietly set obs 100
set seed 10101
gen t=_n
scalar a = 0.48
scalar b = 0.83
scalar c = 0.005
scalar d = 0.014
gen rho = a+(b-a)*runiform()
gen rho2 = (rho)^2
gen s = c+(d-c)*runiform()
gen et = rnormal(0, s)
quietly gen yt= et in 1/2
quietly replace yt=rho*yt[_n-1]-rho*yt[_n-2]+et in 3/L // rho = psi_1
tsset t
ac yt // indicates MA(3)
pac yt // indicates AR(1)