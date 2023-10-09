clear all

import excel "C:\Users\nafew\OneDrive - East West University\Academic\1. Fall-2022\ECO 489 - Forecasting\Assignmnets\1\files\excel\EHL - Daily.xlsx", sheet("stata-data") firstrow

tsset Date, daily

gen l1_return = (Return[_n-1])

gen l2_return = (Return[_n-2])

gen l3_return = (Return[_n-3])

gen l4_return = (Return[_n-4])



reg Return l1_return
outreg2 using As-1.doc, replace ctitle(AR 1)

reg Return l1_return l2_return
outreg2 using As-1.doc, append ctitle(AR 2)

reg Return l1_return l2_return l3_return l4_return
outreg2 using As-1.doc, append ctitle(AR 4)