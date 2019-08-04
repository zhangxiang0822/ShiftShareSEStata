version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

log using "output/logfiles/preclean.log", replace

use "data/ADH_emp_share.dta", clear
renvars, subst(var emp_share)
save "temp/ADH_emp_share.dta", replace

use "data/ADHdata.dta", clear
merge 1:1 czone year using "temp/ADH_emp_share.dta", assert(3) nogen
sort year czone

save "data/ADH_derived.dta", replace

cap log close

