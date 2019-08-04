version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

log using "output/logfiles/AKM_OLS.log", replace

prog main
/*
	**** AKM OLS examples
	** Compute AKM standard error and confidence interval
	use "data/ADH_derived.dta", clear
	
	AKM_OLS, dependant_var(d_sh_empl) shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95) akmtype(1) 
			 
	** Compute AKM0 standard error and confidence interval
	use "data/ADH_derived.dta", clear
	AKM_OLS, dependant_var(d_sh_empl) shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
			 share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95) akmtype(0) 
			 
	** Compute AKM standard error and confidence interval
	use "data/ADH_derived.dta", clear
	
	AKM_OLS, dependant_var(d_sh_empl_mfg) shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95) akmtype(1) 
			 
	** Compute AKM0 standard error and confidence interval
	use "data/ADH_derived.dta", clear
	AKM_OLS, dependant_var(d_sh_empl_mfg) shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
			 share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95) akmtype(0) 
	
	**** AKM IV examples
	** Compute AKM standard error and confidence interval
	use "data/ADH_derived.dta", clear
	
	AKM_IV,  dependant_var(d_sh_empl_mfg) endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95) akmtype(1) 
	*/
	** Compute AKM0 standard error and confidence interval
	use "data/ADH_derived.dta", clear
	
	AKM_IV,  dependant_var(d_sh_empl_mfg) endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95) akmtype(1) 
	/*
	use "data/ADH_derived.dta", clear
	AKM_IV,  dependant_var(d_sh_empl_mfg) endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95) akmtype(0) beta0(0) 
			 */
end

main
