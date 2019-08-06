version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

set matsize 10000
use "data/ADH_derived.dta", clear
preserve
	
use "data/sector_derived.dta", clear
tostring sector, replace

levelsof sec_3d
local sector_3d_list = r(levels)
display `sector_3d_list'
local count = 1

foreach sector_3d in `sector_3d_list' {
	display `sector'
	
	* local sector`count' = ""
	local haha = ""
	qui levelsof sector if sec_3d == `sector_3d'
	
	local sec_list`count' = r(levels)
	local count = `count' + 1
}

restore
local count = `count' - 1

display `sec_list1'
local small_count = 1

foreach sector in `sec_list1'{
	mkmat emp_share`sector', matrix(temp`small_count')
	local small_count = `small_count' + 1
} 
local small_count = `small_count' - 1

matrix small_share = temp1
forvalues i = 2(1)`small_count' {
	matrix small_share = small_share, temp`i'
}
mat list small_share
