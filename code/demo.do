version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

set matsize 10000

use "data/ADH_derived.dta", clear
keep emp_share* 
mkmat emp_share*, matrix(W)

use "data/sector_derived.dta", clear

sort sec_3d sec
tostring sec_3d, replace
	
qui describe
local num_sector = r(N)
			
local current_cluster = ""
			
forvalues i = 1(1)`num_sector' {
	local sec_3d = sec_3d[`i']
	display `i'
	if ("`sec_3d'" != "`current_cluster'") {
		local current_cluster = "`sec_3d'"
		
		matrix share_matrix = W[1..., `i']
	}
	else {
		matrix temp_share_matrix = W[1..., `i']
		matrix share_matrix = share_matrix, temp_share_matrix

	}
}

		
		matrix list share_matrix
