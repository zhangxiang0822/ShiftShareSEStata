version 14.0
clear all
set more off
set scheme s1color
capture log close

cd "D:/Dropbox/RA_Princeton/RA_Morales/BartikSEStata"

log using "output/logfiles/AKM_OLS.log", replace

prog main

	** Compute AKM standard error and confidence interval
	use "data/ADH_derived.dta", clear
	
	AKM_nocluster, dependant_var(d_sh_empl) shiftshare_var(d_tradeotch_pw_lag) ///
				   control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
				   share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95)
	
	** Compute AKM0 standard error and confidence interval
	use "data/ADH_derived.dta", clear
	AKM0_nocluster, dependant_var(d_sh_empl) shiftshare_var(d_tradeotch_pw_lag) ///
				   control_varlist(t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource) ///
				   share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.95)
end

prog AKM_nocluster
	syntax, dependant_var(str) shiftshare_var(str) share_varlist(str) alpha(str) [control_varlist(str) weight_var(str)]
	        
	set matsize 10000
	
	** Generate constant term
	qui gen constant = 1
	
	** Reweight variables according to the weight if weights are given
	if ("`weight_var'" ~= "") {
		foreach var in `dependant_var' `shiftshare_var' constant  {
			qui replace `var' = `var' * sqrt(`weight_var')
		}
		
		foreach var of varlist `share_varlist' {
			qui replace `var' = `var' * sqrt(`weight_var')
		}
		
		** If there are control vars
		if ("`control_varlist'" ~= "") {
			foreach var in `control_varlist' {
				qui replace `var' = `var' * sqrt(`weight_var')
			}
		}
	}
	
	/*
	* Drop colinear share variables
	_rmcoll emp_share*, force
	local share_emp_vars `r(varlist)'
	keep `r(varlist)'  `dependant_var' `shiftshare_var' `control_varlist' constant czone year `weight_var'
	*/
	
	** Generate Matrix of Regressors, shares, and outcome variable
	if ("`control_varlist'" ~= "") {
		mkmat `shiftshare_var' `control_varlist' constant, matrix(Mn)   //Matrix of regressors
	}
	else {
		mkmat `shiftshare_var' constant, matrix(Mn)     //Matrix of regressors
	}
	mkmat emp_share*, matrix(ln)						//Matrix of Shares
	mkmat `dependant_var', matrix(tildeYn) 				//Dependent Variable  
	
	** OLS Estimates
	mat hat_theta = inv(Mn'*Mn)*(Mn' * tildeYn)
	mat e = tildeYn - Mn * hat_theta
	mat hat_beta = hat_theta[1, .]
	local coef = hat_theta[1,1]
    
	** Auxiliary variables
	mkmat `control_varlist' constant, matrix(tildeZn)
	mkmat `shiftshare_var', matrix(tildeXn)
	local dim = rowsof(tildeXn)
	
	mat I_mat = I(`dim')
	mat Ydd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeYn
	mat Xdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeXn
	mat Xddd = inv(ln'*ln) * (ln' * Xdd)
	
	** Compute SE
	mat R_raw = (e' * ln)'
	svmat R_raw, names(R_raw)
	svmat Xddd, names(Xddd)
	qui drop if mi(R_raw)
	
	gen R_raw_sq = R_raw^2
	gen Xddd_sq = Xddd^2

	mkmat R_raw_sq, matrix(R_sq)
	mkmat Xddd_sq, matrix(Xddd_sq)
	
	mat LambdaAKM = R_sq' * Xddd_sq

	mat variance = inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
	local SE_AKM = sqrt(variance[1,1])
	local critical_value = invnormal(0.5 + `alpha'/2)
	local CI_low = `coef' - `critical_value' * `SE_AKM'
	local CI_upp = `coef' + `critical_value' * `SE_AKM'
	
	local tstat = `coef' / `SE_AKM'
	local p = 2*(1 - normal(abs(`tstat')))
	
	if ("`weight_var'" ~= "") {
		qui reg `dependant_var' `shiftshare_var' `control_varlist' [aw = `weight_var']
		
		local SE_homo = _se[`shiftshare_var']
		local t_homo  = _b[`shiftshare_var']/_se[`shiftshare_var']
		local p_homo  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_homo = _b[`shiftshare_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`shiftshare_var'] + `critical_value' * `SE_homo'
		
		qui reg `dependant_var' `shiftshare_var' `control_varlist' [aw = `weight_var'], r
		local SE_r = _se[`shiftshare_var']
		local t_r  = _b[`shiftshare_var']/_se[`shiftshare_var']
		local p_r  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_r = _b[`shiftshare_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`shiftshare_var'] + `critical_value' * `SE_r'
		
		* reg `dependant_var' `shiftshare_var' `control_varlist' [aw = `weight_var'], cluster(state)
	}
	else {
		qui reg `dependant_var' `shiftshare_var' `control_varlist'
		local SE_homo = _se[`shiftshare_var']
		local t_homo  = _b[`shiftshare_var']/_se[`shiftshare_var']
		local p_homo  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_homo = _b[`shiftshare_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`shiftshare_var'] + `critical_value' * `SE_homo'
		
		qui reg `dependant_var' `shiftshare_var' `control_varlist', r
		local SE_r = _se[`shiftshare_var']
		local t_r  = _b[`shiftshare_var']/_se[`shiftshare_var']
		local p_r  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_r = _b[`shiftshare_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`shiftshare_var'] + `critical_value' * `SE_r'
		
		* reg `dependant_var' `shiftshare_var' `control_varlist', cluster(state)
	}
	
	** Output Results
	display "The estimated coefficient is " %5.4f `coef'
	display "Inference"
	display "               Std. Error   p-value   Lower CI   Upper CI"
	display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
	display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
	display "AKM               " %5.4f `SE_AKM'   "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp' 
end


prog AKM0_nocluster
	syntax, dependant_var(str) shiftshare_var(str) share_varlist(str) alpha(str) [control_varlist(str) weight_var(str)]
	
	set matsize 10000

	** Generate constant term
	qui gen constant = 1
	
	** Reweight variables according to the weight if weights are given
	if ("`weight_var'" ~= "") {
		foreach var in `dependant_var' `shiftshare_var' constant  {
			qui replace `var' = `var' * sqrt(`weight_var')
		}
		
		foreach var of varlist `share_varlist' {
			qui replace `var' = `var' * sqrt(`weight_var')
		}
		
		** If there are control vars
		if ("`control_varlist'" ~= "") {
			foreach var in `control_varlist' {
				qui replace `var' = `var' * sqrt(`weight_var')
			}
		}
	}
	
	** Generate Matrix of Regressors, shares, and outcome variable
	if ("`control_varlist'" ~= "") {
		mkmat `shiftshare_var' `control_varlist' constant, matrix(Mn)   //Matrix of regressors
	}
	else {
		mkmat `shiftshare_var' constant, matrix(Mn)     //Matrix of regressors
	}
	mkmat emp_share*, matrix(ln)						//Matrix of Shares
	mkmat `dependant_var', matrix(tildeYn) 				//Dependent Variable  
	
	** OLS Estimates
	mat hat_theta = inv(Mn'*Mn)*(Mn' * tildeYn)
	mat e = tildeYn - Mn * hat_theta
	mat hat_beta = hat_theta[1, .]
	local coef = hat_theta[1,1]
    
	** Auxiliary variables
	mkmat `control_varlist' constant, matrix(tildeZn)
	mkmat `shiftshare_var', matrix(tildeXn)
	local dim = rowsof(tildeXn)
	
	mat I_mat = I(`dim')
	mat Ydd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeYn
	mat Xdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeXn
	mat Xddd = inv(ln'*ln) * (ln' * Xdd)
	
	** Compute SE
	local beta0 = 0
	mat e_null = Ydd - Xdd * `beta0'
	
	mat R_raw = (e_null' * ln)'
	svmat R_raw, names(R_raw)
	svmat Xddd, names(Xddd)
	qui drop if mi(R_raw)
	
	gen R_raw_sq = R_raw^2
	gen Xddd_sq = Xddd^2

	mkmat R_raw_sq, matrix(R_sq)
	mkmat Xddd_sq, matrix(Xddd_sq)
	
	mat LambdaAKM = R_sq' * Xddd_sq

    * Variance matrix
	mat variance = inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
	local SE_AKMnull_n = sqrt(variance[1,1])
	
	* Compute Confidence INterval
	local critical_value = invnormal(1-0.05/2)
    local critical2 = `critical_value'^2
    mat RY = Xdd' * Ydd
    mat RX = Xdd' * Xdd
	
	mat lnY =  (Ydd' * ln)'
    mat lnX =  (Xdd' * ln)'
	
	svmat lnY, names(lnY)
	svmat lnX, names(lnX)
	gen lnY_lnX = lnY * lnX
	gen lnY_lnY = lnY1 * lnY1
	gen lnX_lnX = lnX1 * lnX1
	
    mkmat lnY_lnX, matrix(lnY_lnX)
    mkmat lnY_lnY, matrix(lnY_lnY)
    mkmat lnX_lnX, matrix(lnX_lnX)
	
    mat SXY = lnY_lnX' * Xddd_sq
    mat SXX = lnX_lnX' * Xddd_sq
    mat SYY = lnY_lnY' * Xddd_sq
	
	mat Q = (RX * RX)/`critical2' - SXX
    mat Delta = (RY * RX - `critical2' * SXY) * (RY * RX - `critical2' * SXY) - (RX * RX - `critical2' * SXX) * (RY * RY - `critical2' * SYY)
	
	local Q = Q[1,1]
	local Delta = Delta[1,1]
	if `Q' > 0 {
		mat CIl = ((RY * RX - `critical2' * SXY) - `Delta'^(1/2) ) * inv(RX * RX - `critical2' * SXX)
		mat CIu = ((RY * RX - `critical2' * SXY) + `Delta'^(1/2) ) * inv(RX * RX - `critical2' * SXX)
		local CI_low = CIl[1,1]
		local CI_upp = CIu[1,1]
		local CIType = 1
	} 
	else {
		if `delta' > 0 {
			mat CIl = ((RY*RX - `critical2' * SXY) + `Delta'^(1/2) ) * inv(RX * RX - `critical2' * SXX)
			mat CIu = ((RY*RX - `critical2' * SXY) - `Delta'^(1/2) ) * inv(RX * RX - `critical2' * SXX)
			local CI_low = CIl[1,1]
			local CI_upp = CIu[1,1]
			local CIType = 2
		} 
		else {
			local CI_low = -10000000
			local CI_upp = 10000000
			local CIType = 3
		}
	}    
    
	local SE_AKM0 = (`CI_upp' - `CI_low')/(2 * `critical_value')
	local tstat = (`coef' - `beta0') / `SE_AKMnull_n'
	local p = 2*(1 - normal(abs(`tstat')))
	
	if ("`weight_var'" ~= "") {
		qui reg `dependant_var' `shiftshare_var' `control_varlist' constant [aw = `weight_var'], noconstant
		
		local SE_homo = _se[`shiftshare_var']
		local t_homo  = _b[`shiftshare_var']/_se[`shiftshare_var'] 
		local p_homo  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_homo = _b[`shiftshare_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`shiftshare_var'] + `critical_value' * `SE_homo'
		
		qui reg `dependant_var' `shiftshare_var' `control_varlist' constant [aw = `weight_var'], r noconstant
		local SE_r = _se[`shiftshare_var']
		local t_r  = _b[`shiftshare_var']/_se[`shiftshare_var']
		local p_r  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_r = _b[`shiftshare_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`shiftshare_var'] + `critical_value' * `SE_r'
		
		* reg `dependant_var' `shiftshare_var' `control_varlist' [aw = `weight_var'], cluster(state)
	}
	else {
		qui reg `dependant_var' `shiftshare_var' `control_varlist' constant, noconstant
		local SE_homo = _se[`shiftshare_var']
		local t_homo  = _b[`shiftshare_var']/_se[`shiftshare_var']
		local p_homo  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_homo = _b[`shiftshare_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`shiftshare_var'] + `critical_value' * `SE_homo'
		
		qui reg `dependant_var' `shiftshare_var' `control_varlist'  constant, r noconstant
		local SE_r = _se[`shiftshare_var']
		local t_r  = _b[`shiftshare_var']/_se[`shiftshare_var']
		local p_r  = 2*ttail(e(df_r),abs(`t_homo'))
		local CI_low_r = _b[`shiftshare_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`shiftshare_var'] + `critical_value' * `SE_r'
		
		* reg `dependant_var' `shiftshare_var' `control_varlist', cluster(state)
	}
	
	** Output Results
	display "The estimated coefficient is " %5.4f `coef'
	display "Inference"
	display "               Std. Error   p-value   Lower CI   Upper CI"
	display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
	display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
	display "AKM0              " %5.4f `SE_AKM0'  "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp' 
end

main

cap log close
