program define ivreg_ss, eclass
	version 14.0
	syntax varlist(numeric min=1 max=1), endogenous_var(varlist numeric min=1 max=1) ///
		   shiftshare_iv(varlist numeric min=1 max=1) share_varlist(varlist numeric min=1) ///
		   [control_varlist(varlist numeric min=1) weight_var(varlist numeric min=1 max=1) ///
		    akmtype(str) beta0(real 0.0) alpha(str) path_cluster(str) cluster_var(str) firststage(integer 0)]
	
	set more off
	set matsize 10000
	
	local dependant_var `varlist'
	
	if "`beta0'" == "" {
		local beta0 = 0
	}
	
	if "`alpha'" == "" {
		local alpha = 0.05
	}
	
	if "`akmtype'" == "" {
		local akmtype = 1
	}
	
	
	** Show first-stage results
	if `firststage' != 0 {
		display ""
		display "Below we show First-stage results"
		if "`control_varlist'" ~= "" {
			if "`weight_var'" ~= "" {
				if "`cluster_var'" ~= "" {
					display "We're sorry that we don't support reporting first-stage results with clustering. To obtain the correct clustered AKM standard errors, please use reg_ss command to run first-stage regression manually."
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
						 alpha(`alpha') control_varlist(`control_varlist') weight_var(`weight_var') akmtype(`akmtype')
				}
				else {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
						 alpha(`alpha') control_varlist(`control_varlist') weight_var(`weight_var') akmtype(`akmtype')
				}
			}
			else {
				if "`cluster_var'" ~= "" {
					display "We're sorry that we don't support reporting first-stage results with clustering. To obtain the correct clustered AKM standard errors, please use reg_ss command to run first-stage regression manually."
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') control_varlist(`control_varlist') akmtype(`akmtype')
				}
				else {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') control_varlist(`control_varlist') akmtype(`akmtype')
				}
			}
		}
		else {
			if "`weight_var'" ~= "" {
				if "`cluster_var'" ~= "" {
					display "We're sorry that we don't support reporting first-stage results with clustering. To obtain the correct clustered AKM standard errors, please use reg_ss command to run first-stage regression manually."
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') weight_var(`weight_var') akmtype(`akmtype') 
				}
				else {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') weight_var(`weight_var') akmtype(`akmtype')
				}
			}
			else {
				if "`cluster_var'" ~= "" {
					display "We're sorry that we don't support reporting first-stage results with clustering. To obtain the correct clustered AKM standard errors, please use reg_ss command to run first-stage regression manually."
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') akmtype(`akmtype')
				}
				else {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') akmtype(`akmtype')
				}
			}
		}
		
		local b_firststage = `e(b)'
		local se_firststage = `e(se)'
		local CI_upp_firststage = `e(CI_upp)'
		local CI_low_firststage = `e(CI_low)'
		local p_firststage = `e(p)'
		local tstat_firststage = `e(tstat)'
	}

	** Some locals
	local critical_value = invnormal(1- `alpha'/2)
	
	if ("`weight_var'" ~= "") {
		qui ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var']
		local SE_homo = _se[`endogenous_var']
		local z_homo  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_homo  = normal(`z_homo')
		local CI_low_homo = _b[`endogenous_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`endogenous_var'] + `critical_value' * `SE_homo'
		
		qui ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], r 
		local SE_r = _se[`endogenous_var']
		local z_r  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_r  = normal(`z_r')
		local CI_low_r = _b[`endogenous_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`endogenous_var'] + `critical_value' * `SE_r'
		
		* ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], cluster(state)
	}
	else {
		qui ivregress 2sls `dependant_var' `control_varlist'  (`endogenous_var' = `shiftshare_iv') 
		local SE_homo = _se[`endogenous_var']
		local z_homo  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_homo  = normal(`z_homo')
		local CI_low_homo = _b[`endogenous_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`endogenous_var'] + `critical_value' * `SE_homo'
		
		qui ivregress 2sls `dependant_var' `control_varlist'  (`endogenous_var' = `shiftshare_iv'), vce(robust)
		local SE_r = _se[`endogenous_var']
		
		local z_r  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_r  = normal(`z_r')
		local CI_low_r = _b[`endogenous_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`endogenous_var'] + `critical_value' * `SE_r'
		
		* ivregress 2sls `dependant_var' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], cluster(state)
	}
	
	** Generate constant term
	qui gen constant = 1
	
	** Reweight variables according to the weight if weights are given
	if ("`weight_var'" ~= "") {
		foreach var in `dependant_var' `endogenous_var' `shiftshare_iv' constant  {
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
	
	capture _rmcoll `control_varlist' constant, force nocons
	if _rc == 0{
		_rmcoll `control_varlist' constant, force nocons
		local controls `r(varlist)'
	}
	else {
		display "Warning: You must manually generate dummy variables, instead of using i.XXX"
		exit
	}

	_rmcoll `share_varlist', force
	local num_omit_var = `r(k_omitted)'
	
	if `num_omit_var' > 0 {
		display "Error: You have collinear share variables (Share matrix has colinear columns)"
		exit
	}
	
	** Generate Matrix of Regressors, shares, and outcome variable
	if ("`control_varlist'" ~= "") {
		mkmat `shiftshare_iv' `controls', matrix(Mn) //Matrix of IVs
		mkmat `endogenous_var' `controls', matrix(Gn)   //Matrix of regressors
	}
	else {
		mkmat `shiftshare_iv' constant, matrix(Mn)     //Matrix of IVs
		mkmat `endogenous_var' constant, matrix(Gn)    //Matrix of regressors
	}
	mkmat `share_varlist', matrix(ln)						//Matrix of Shares
	mkmat `dependant_var', matrix(tildeYn) 				//Dependent Variable   
	
	local num_counties = rowsof(ln)
	local num_sector1  = colsof(ln)
	
	if (`num_counties' < `num_sector1') {
		display "ERROR: You have more number of sectors than regions"
		exit
	}
	
	** OLS Estimates
	mat P1 = Mn' * Gn
	mat P2 = inv(Mn' * Mn)
	mat P3 = Mn' * tildeYn
	mat hat_theta = inv(P1'*P2*P1) * (P1'*P2*P3)
	
	mat e = tildeYn - Gn * hat_theta
	local hat_beta = hat_theta[1, 1]
	
	** Auxiliary variables
	mkmat `controls', matrix(tildeZn)
	mkmat `shiftshare_iv', matrix(tildeXn)
	mkmat `endogenous_var', matrix(tildeGn)
	
	mat A = tildeZn * inv(tildeZn'*tildeZn)
	mat Xdd = tildeXn - A*(tildeZn'*tildeXn)
	mat Ydd = tildeYn - A*(tildeZn'*tildeYn) 
	mat Gdd = tildeGn - A*(tildeZn'*tildeGn)
	mat Xddd = inv(ln'*ln) * (ln'*Xdd)
	
	** First stage coefficient
	mat hat_thetaFS = inv(Mn'*Mn) * (Mn'*Gdd)
	local hatpi = hat_thetaFS[1,1]
	
	** Estimate AKM standard error
	if `akmtype' == 1 {
		** Compute SE
		if "`path_cluster'" == "" {
			mat R_raw = (e' * ln)'
			svmat R_raw, names(R_raw)
			svmat Xddd, names(Xddd)
			qui drop if mi(R_raw)
			
			gen R_raw_sq = R_raw^2
			gen Xddd_sq = Xddd^2
			
			mkmat R_raw_sq, matrix(R_sq)
			mkmat Xddd_sq, matrix(Xddd_sq)
			
			mat LambdaAKM = R_sq' * Xddd_sq
			
		}
		else {		
			** Get list of share variables by cluster
			use "`path_cluster'", clear
			
			qui sum `cluster_var'
			local num_sector2 = r(N)
			
			if (`num_sector1' != `num_sector2') {
				display "Error: The length of your cluster_var is different from the number of sectors"
				exit
			}
			
			mat e_ln = (e' * ln)'
			svmat e_ln, names(e_ln)
			svmat Xddd, names(Xddd)
			qui keep if ~mi(e_ln)
			sort `cluster_var'
		
			matrix opaccum A = e_ln, group(`cluster_var') opvar(Xddd)
			mat LambdaAKM = A[1,1]
		}
		mat variance = 1/(`hatpi')^2 * inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
		local SE_AKM = sqrt(variance[1,1])
		
		local CI_low = `hat_beta' - `critical_value' * `SE_AKM'
		local CI_upp = `hat_beta' + `critical_value' * `SE_AKM'
		
		local tstat = `hat_beta' / `SE_AKM'
		local p = 2*(1 - normal(abs(`tstat')))
		
		** Output Results
		display " "
		display "Below we show IV regression results"
		display "The estimated coefficient is " %5.4f `hat_beta'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
		display "AKM               " %5.4f `SE_AKM'   "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp'
	}
		
	if `akmtype' == 0 {
		** Compute SE
		mat e_null = Ydd - Gdd * `beta0'
		
		svmat Xddd, names(Xddd)
		qui gen Xddd_sq = Xddd^2
		
		if "`path_cluster'" == "" {
			mat R_raw = (e_null' * ln)'
			svmat R_raw, names(R_raw)
			
			qui drop if mi(R_raw)
			
			gen R_raw_sq = R_raw^2
			mkmat R_raw_sq, matrix(R_sq)
			mkmat Xddd_sq, matrix(Xddd_sq)
			
			mat LambdaAKM = R_sq' * Xddd_sq
		}
		else {
			** Get list of share variables by cluster
			use "`path_cluster'", clear
			
			qui sum `cluster_var'
			local num_sector2 = r(N)
			
			if (`num_sector1' != `num_sector2') {
				display "Error: The length of your cluster_var is different from the number of sectors"
				exit
			}
			
			mat e_ln = (e_null' * ln)'
			mat Xdd_ln = (Gdd' * ln)'
			mat Ydd_ln = (Ydd' * ln)'
			svmat e_ln, names(e_ln)
			svmat Xdd_ln, names(Xdd_ln)
			svmat Ydd_ln, names(Ydd_ln)
			svmat Xddd, names(Xddd)
			
			qui keep if ~mi(e_ln)
			sort `cluster_var'
		
			matrix opaccum A = e_ln, group(`cluster_var') opvar(Xddd)
			mat LambdaAKM = A[1,1]

			matrix opaccum B = Xdd_ln, group(`cluster_var') opvar(Xddd)
			mat SXX = B[1,1]
			
			matrix opaccum B = Ydd_ln, group(`cluster_var') opvar(Xddd)
			mat SYY = B[1,1]
			
			qui levelsof `cluster_var'
			local sector_list = r(levels)
			
			mat SXY = J(1,1,0)
			foreach cat in `sector_list' {
				preserve
				qui keep if `cluster_var' == `cat'
				
				mkmat Xdd_ln, matrix(Xdd_ln1)
				mkmat Ydd_ln, matrix(Ydd_ln1)
				mkmat Xddd, matrix(Xddd1) 
				
				mat SXY = SXY + Xdd_ln1'*Xddd1*Xddd1'*Ydd_ln1
				restore
			}
		}
		
		* Variance matrix
		mat variance = 1/`hatpi'^2 * inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
		local SE_AKMnull_n = sqrt(variance[1,1])
		
		* Compute Confidence INterval
		local critical_value = invnormal(1-0.05/2)
		local critical2 = `critical_value'^2
		mat RY = Xdd' * Ydd
		mat RX = Xdd' * Gdd

		if "`path_cluster'" == "" {
			mat lnY =  (Ydd' * ln)'
			mat lnX =  (Gdd' * ln)'
			
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
		}

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
			if `Delta' > 0 {
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
		
		local tstat = (`hat_beta' - `beta0') / `SE_AKMnull_n'
		local p = 2*(1 - normal(abs(`tstat')))
		
		** Output Results
		display " "
		display "Below we show IV regression results"
		display "The estimated coefficient is " %5.4f `hat_beta'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
		
		if `CIType' == 1 {
			display "AKM0              " %5.4f `SE_AKM0'  "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp'
		}
		else if `CIType' == 2 {
			display "Inference"
			display "               Std. Error   p-value   CI"
			display "AKM0              " %5.4f `SE_AKM0'  "    " %5.4f `p' "    " "=(-Inf, " %5.4f `CI_low' "] + [" %5.4f `CI_upp' ", Inf)"
		}
		else {
			display "Inference"
			display "               Std. Error   p-value   CI"
			display "AKM0              " %5.4f `SE_AKM0'  "    " %5.4f `p' "    " "=(-Inf, Inf)"
		}
	}	

	
	** Save results for display
	ereturn scalar b = `hat_beta'
	if `akmtype' == 0 {
		ereturn local se = `SE_AKM0'
	}
	else {
		ereturn local se = `SE_AKM'
	}

	ereturn local CI_upp = `CI_upp'
	ereturn local CI_low = `CI_low'
	ereturn scalar p = `p'
	ereturn scalar tstat = `tstat'
	
	if `firststage' != 0 {
		ereturn local se_firststage = `se_firststage'
		ereturn local b_firststage  = `b_firststage'
		ereturn local CI_upp_firststage = `CI_upp_firststage'
		ereturn local CI_low_firststage = `CI_low_firststage'
		ereturn local p_firststage = `p_firststage'
		ereturn scalar tstat = `tstat_firststage'
	}
	*/
end
