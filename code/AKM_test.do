version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

log using "output/logfiles/AKM_test.log", replace

/*
This is the test dofile for Stata package ivreg_ss and reg_ss. The corresponding
R version test can be found at:
	- https://github.com/kolesarm/ShiftShareSE/blob/master/tests/testthat/test_se.R
*/

prog main	

	* c1
	use "data/ADH_derived.dta", clear
	reg_ss d_tradeusch_pw, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)

	local p_c1_akm = `e(p)'
	local ci_c1_low_akm = `e(CI_low)'
	local ci_c1_upp_akm = `e(CI_upp)'

	use "data/ADH_derived.dta", clear
	reg_ss d_tradeusch_pw, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	local c1_coef = `e(b)'
	local p_c1_akm0 = `e(p)'
	local ci_c1_low_akm0 = `e(CI_low)'
	local ci_c1_upp_akm0 = `e(CI_upp)'
	
	* c3
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	local p_c3_akm = `e(p)'
	local ci_c3_low_akm = `e(CI_low)'
	local ci_c3_upp_akm = `e(CI_upp)'
	
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	local c3_coef = `e(b)'
	local p_c3 = `e(p)'
	local p_c3_akm0 = `e(p)'
	local ci_c3_low_akm0 = `e(CI_low)'
	local ci_c3_upp_akm0 = `e(CI_upp)'

	* c5
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)

	local p_c5_akm = `e(p)'
	local ci_c5_low_akm = `e(CI_low)'
	local ci_c5_upp_akm = `e(CI_upp)'
	
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	local c5_coef = `e(b)'
	local p_c5 = `e(p)'
	local p_c5_akm0 = `e(p)'
	local ci_c5_low_akm0 = `e(CI_low)'
	local ci_c5_upp_akm0 = `e(CI_upp)'
		 
	* b1
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
		 
	* b3
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	local p_b3 = `e(p)'
	
	* b5
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	
	local p_b5 = `e(p)'
	
	* a1
	use "data/ADH_derived.dta", clear
	reg_ss d_tradeusch_pw, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) 
	
	local p_a1_akm = `e(p)'
	local ci_a1_low_akm = `e(CI_low)'
	local ci_a1_upp_akm = `e(CI_upp)'
	
	use "data/ADH_derived.dta", clear
	reg_ss d_tradeusch_pw, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) 
	
	local a1_coef = `e(b)'
	local p_a1_akm0 = `e(p)'
	local ci_a1_low_akm0 = `e(CI_low)'
	local ci_a1_upp_akm0 = `e(CI_upp)'
	
	* a3 
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) 
	
	local p_a3_akm = `e(p)'
	local ci_a3_low_akm = `e(CI_low)'
	local ci_a3_upp_akm = `e(CI_upp)'
	
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) 
			 
	local p_a3 = `e(p)'
	local p_a3_akm0 = `e(p)'
	local ci_a3_low_akm0 = `e(CI_low)'
	local ci_a3_upp_akm0 = `e(CI_upp)'
	
	* a5
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) 
	
	local p_a5_akm = `e(p)'
	local ci_a5_low_akm = `e(CI_low)'
	local ci_a5_upp_akm = `e(CI_upp)'
	
	
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) 
	
	local p_a5 = `e(p)'
	local p_a5_akm0 = `e(p)'
	local ci_a5_low_akm0 = `e(CI_low)'
	local ci_a5_upp_akm0 = `e(CI_upp)'
	
	** Test
	** (1) b3 p-value = b5 p-value
	display ""
	display "b3 and b5 AKM0 p-value equal"
	display `p_b3'
	display `p_b5'
	
	** (2) c3 p-value = c5 p-value
	display ""
	display "c3 and c5 AKM0 p-value equal"
	display `p_c3'
	display `p_c5'
	
	** (3) a3 p-value = a5 p-value
	display ""
	display "a3 and a5 AKM0 p-value equal"
	display `p_a3'
	display `p_a5'
	
	** (4) a3 and b3 results identical
	** (5) b5 and a5 results identical
	** (6) Check c1-c3 coefficient
	display ""
	display "Check c1-c3 coefficient"
	display "Expected value: 0.631040938185757, -0.488568717122902, -0.774226658776735"
	display "Realized Value: `c1_coef', `c3_coef', `c5_coef' "
	
	** (7) c1 coefficient equals to a1 coefficient
	display ""
	display "c1 coefficient equals to a1 coefficient"
	display "`c1_coef' = `a1_coef'"
	
	** (8) c1 and a1 p-value, CI-low, and CI-upp 
	display ""
	display "c1 and a1 p-value, CI-low, and CI-upp"
	display "Expected value: 0, 0.001282890950014, 0, 0.000056519706974"
    display "Realized value: `p_c1_akm', `p_c1_akm0', `p_a1_akm', `p_a1_akm0'"
	
	display ""
	display "Expected value: 0.527240172127609, 0.537570958010291, 0.495188795569771, 0.522017706924849"
    display "Realized value: `ci_c1_low_akm', `ci_c1_low_akm0', `ci_a1_low_akm', `ci_a1_low_akm0'"
	
	display ""
	display "Expected value: 0.734841704243901, 0.838282656411893, 0.766893080801748, 0.889797428206052"
    display "Realized value: `ci_c1_upp_akm', `ci_c1_upp_akm0', `ci_a1_upp_akm', `ci_a1_upp_akm0'"
	
	** (9) c3 and a3 p-value, CI-low, and CI-upp 
	display ""
	display "c3 and a3 p-value, CI-low, and CI-upp"
	display "Expected value: 0.002924641237718, 0.000421803253836, 0.000000150063370, 0.000090458689727"
    display "Realized value: `p_c3_akm', `p_c3_akm0', `p_a3_akm', `p_a3_akm0'"
	
	display ""
	display "Expected value: -0.810383925354586, -1.236885319089530, -0.516754284756085, -0.629980318872062"
    display "Realized value: `ci_c3_low_akm', `ci_c3_low_akm0', `ci_a3_low_akm', `ci_a3_low_akm0'"
	
	display ""
	display "Expected value: -0.166753508891236, -0.239754057137763, -0.235900929361803, -0.257494610447470"
    display "Realized value: `ci_c3_upp_akm', `ci_c3_upp_akm0', `ci_a3_upp_akm', `ci_a3_upp_akm0'"
	
	** (10) c5 and a5 p-value, CI-low, and CI-upp 
	display ""
	display "c5 and a5 p-value, CI-low, and CI-upp"
	display "Expected value: 0.001277718157450, 0.000421803253837, 0.000000051567025, 0.000090458689726"
    display "Realized value: `p_c5_akm', `p_c5_akm0', `p_a5_akm', `p_a5_akm0'"
	
	display ""
	display "Expected value: -1.245349169792894, -1.690324047107167, -0.810991466572728, -0.891427381124374"
    display "Realized value: ``ci_c5_low_akm', `ci_c5_low_akm0', `ci_a5_low_akm', `ci_a5_low_akm0'"
	
	display ""
	display "Expected value: -0.303104147760458, -0.389313219498141, -0.381728638531380, -0.391771441695858"
    display "Realized value: ``ci_c5_upp_akm', `ci_c5_upp_akm0', `ci_a5_upp_akm', `ci_a5_upp_akm0'"

	** Test_that: Weak IV
	** iv0	
	use "data/ADH_derived.dta", clear
	gen division =  2*reg_midatl + 3*reg_encen + 4*reg_wncen + 5*reg_satl+ 6*reg_escen + 7*reg_wscen + 8*reg_mount + 9 * reg_pacif
	keep if division < 8
	
	forvalues i = 1(1)7 {
		gen division`i' = (division == `i')
	}

	local control_varlist t2 l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource division1-division7
	
	ivreg_ss d_sh_empl, endogenous_var(d_tradeusch_pw) ///
						shiftshare_iv(d_tradeotch_pw_lag) ///
						control_varlist(`control_varlist') ///
						share_varlist(emp_share1-emp_share770) ///
						alpha(0.05) akmtype(0) 
						
	local se_iv0 = `e(se)'
	
	** iv1
	use "data/ADH_derived.dta", clear
	gen division =  2*reg_midatl + 3*reg_encen + 4*reg_wncen + 5*reg_satl+ 6*reg_escen + 7*reg_wscen + 8*reg_mount + 9 * reg_pacif
	keep if division < 7
	
	forvalues i = 1(1)6 {
		gen division`i' = (division == `i')
	}
	local control_varlist t2 l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource division1-division6
	
	ivreg_ss d_sh_empl, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(`control_varlist') ///
		     share_varlist(emp_share1-emp_share770) alpha(0.05) akmtype(0) 

	local se_iv0 = `e(se)'

	** iv2
	use "data/ADH_derived.dta", clear
	gen division =  2*reg_midatl + 3*reg_encen + 4*reg_wncen + 5*reg_satl+ 6*reg_escen + 7*reg_wscen + 8*reg_mount + 9 * reg_pacif
	keep if division < 6
	
	forvalues i = 1(1)5 {
		gen division`i' = (division == `i')
	}
	local control_varlist t2 l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource division1-division5
	
	ivreg_ss d_sh_empl, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(`control_varlist') ///
		     share_varlist(emp_share1-emp_share770) alpha(0.05) akmtype(0) 

	* r0
	use "data/ADH_derived.dta", clear
	gen division =  2*reg_midatl + 3*reg_encen + 4*reg_wncen + 5*reg_satl+ 6*reg_escen + 7*reg_wscen + 8*reg_mount + 9 * reg_pacif
	keep if division >4
	
	forvalues i = 5(1)9 {
		gen division`i' = (division == `i')
	}
	local control_varlist t2 l_shind_manuf_cbp reg_encen reg_midatl reg_satl reg_wncen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource division5-division9
	
	reg_ss d_tradeusch_pw, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(`control_varlist') ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) 

	* r1
	use "data/ADH_derived.dta", clear
	gen division =  2*reg_midatl + 3*reg_encen + 4*reg_wncen + 5*reg_satl+ 6*reg_escen + 7*reg_wscen + 8*reg_mount + 9 * reg_pacif
	keep if division > 4
	
	forvalues i = 5(1)9 {
		gen division`i' = (division == `i')
	}
	local control_varlist t2 l_shind_manuf_cbp reg_encen reg_midatl reg_satl reg_wncen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource division5-division9

	reg_ss d_tradeusch_pw, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(`control_varlist') ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.045) akmtype(0) 
	
	* r2
	use "data/ADH_derived.dta", clear
	gen division =  2*reg_midatl + 3*reg_encen + 4*reg_wncen + 5*reg_satl+ 6*reg_escen + 7*reg_wscen + 8*reg_mount + 9 * reg_pacif
	keep if division <6
	
	forvalues i = 1(1)5 {
		gen division`i' = (division == `i')
	}
	
	local control_varlist t2 l_shind_manuf_cbp l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource division1-division5

	reg_ss d_tradeusch_pw, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(`control_varlist') ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0)
			 
end

cap log close

main
