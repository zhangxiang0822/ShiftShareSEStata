version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

log using "output/logfiles/AKM_OLS.log", replace

prog main
	use "data/ADHdata.dta", clear
	merge 1:1 czone year using "data/ADH_emp_share.dta", assert(3) nogen
	sort year czone
	
	AKM_nocluster
	* AKM0_nocluster
end

prog AKM_nocluster
	local control_var "t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource"

	* (AKM)
	set matsize 10000
	

end

main

cap log close
