version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

log using "output/logfiles/AKM.log", replace

prog main
	**** AKM OLS examples	
	* (1) AKM, no cluster
	timer on 1
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05)
	timer off 1
	timer list 1
	
	* (2) AKM0, no cluster
	timer on 2
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
			 share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) 
	timer off 2
	timer list 2

	* (3) AKM, cluster
	timer on 3
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	timer off 3
	timer list 3
	
	* (4) AKM0, cluster
	timer on 4
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl_mfg, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
			 share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d)
	timer off 4
	timer list 4
	
	**** AKM IV examples
	* (1) AKM, no cluster
	timer on 5
	use "data/ADH_derived.dta", clear
	
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) firststage(1)
	timer off 5
	timer list 5
	
	* (2) AKM0, no cluster
	timer on 6
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) beta0(0) firststage(1)
	timer off 6
	timer list 6
	
	* (3) AKM, cluster
	timer on 7
	use "data/ADH_derived.dta", clear
	
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) path_cluster("data/sector_derived.dta") cluster_var(sec_3d) firststage(1)
	timer off 7
	timer list 7

	* (4) AKM0, no cluster
	timer on 8
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) beta0(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d) firststage(1)
	timer off 8
	timer list 8
	
	** Use d_sh_empl as outcome
	* (1) OLS AKM cluster
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) path_cluster("data/sector_derived.dta") cluster_var(sec_3d) 
	
	* (2） OLS AKM0 cluster
	use "data/ADH_derived.dta", clear
	reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d) 
			 
	* (3) IV AKM cluster
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) beta0(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d) firststage(1)

	* (4) IV AKM0 cluster
	use "data/ADH_derived.dta", clear
	ivreg_ss d_sh_empl, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) ///
			 control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
		     share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) beta0(0) path_cluster("data/sector_derived.dta") cluster_var(sec_3d) firststage(1)
end

cap log close

main
