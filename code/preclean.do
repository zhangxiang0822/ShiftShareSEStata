version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

log using "output/logfiles/preclean.log", replace

use "data/sector_shock.dta", clear
renvars, subst(var sic) 
expand 1444
save "temp/sector_shock.dta", replace

/*
use "data/ADH_emp_share.dta", clear
renvars, subst(var emp_share)
merge 1:1 _n using "temp/sector_shock.dta", assert(3) nogen

** Generate X (shift-share variable)
forvalues i = 1(1)775 {
	gen temp_weights`i' = emp_share`i' * sic`i'
}
egen X = rowtotal(temp_weights*)
drop temp_weights* emp_share* sic*

merge 1:1 czone year using "data/ADHdata.dta", assert(3) nogen
save "data/ADH_derived.dta", replace

cap log close

