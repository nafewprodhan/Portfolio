**********************************************
** ARIMA estimation of Finland real GDP
clear all
set more off

freduse CLVMNACSCAB1GQFI

** set time variable
generate time = q(1990q1) + _n -1
format time %tq
tsset time
rename CLVMNACSCAB1GQFI rgdp
gen lrgdp = log(rgdp)
tsline rgdp
tsline d.rgdp

** Step 1: Identification
ac d.lrgdp if time < tq(2020q1) // MA(3)
pac d.lrgdp if time < tq(2020q1) // AR (3)


** Step 2: Estimation arima(p,d,q)
arima d.lrgdp if time < tq(2020q1), arima(3,0,3)


// arima d.lrgdp, arima(1,0,1)
//
// arima d.lrgdp, arima(2,0,1)

*** Step 3: Diagnostic checking
predict resid, residual
wntestq resid
wntestb resid

** Step 4: Forecasting
tsappend, add(2)
arima d.lrgdp, arima(3,0,3), if tin(,2019q3)
arima d.lrgdp, arima(3,0,3)

// drop lrgdp_f
predict lrgdp_fc, y dynamic(q(2022q4))
predict lrgdp_f, y dynamic(q(2019q4))

** Automatic ARIMA Selection
arimasel d.lrgdp if time < tq(2020q1), ar(4) ma(4)

exit
**********************************************