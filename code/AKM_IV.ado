program define AKM_IV, eclass
	syntax, dependant_var(str) endogenous_var(str) shiftshare_iv(str) share_varlist(str) alpha(str) ///
			akmtype(str) [control_varlist(str) weight_var(str) beta0(str)]
	
	set more off
	set matsize 10000
	
	** IV results withoud AKM adjustment
	local critical_value = invnormal(0.5 + `alpha'/2)
	if ("`weight_var'" ~= "") {
		qui ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var']
		local SE_homo = _se[`endogenous_var']
		local t_homo  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_homo  = 2 * ttail(e(df_r),abs(`t_homo'))
		local CI_low_homo = _b[`endogenous_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`endogenous_var'] + `critical_value' * `SE_homo'
		
		qui ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], r 
		local SE_r = _se[`endogenous_var']
		local t_r  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_r  = 2 * ttail(e(df_r),abs(`t_r'))
		local CI_low_r = _b[`endogenous_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`endogenous_var'] + `critical_value' * `SE_r'
		
		* ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], cluster(state)
	}
	else {
		qui ivregress 2sls `dependant_var' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var']
		local SE_homo = _se[`endogenous_var']
		local t_homo  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_homo  = 2 * ttail(e(df_r),abs(`t_homo'))
		local CI_low_homo = _b[`endogenous_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`endogenous_var'] + `critical_value' * `SE_homo'
		
		qui ivregress 2sls `dependant_var' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], r 
		local SE_r = _se[`endogenous_var']
		local t_r  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_r  = 2 * ttail(e(df_r),abs(`t_r'))
		local CI_low_r = _b[`endogenous_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`endogenous_var'] + `critical_value' * `SE_r'
		
		* ivregress 2sls `dependant_var' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], cluster(state)
	}
		
	if "`beta0'" == "" {
		local beta0 = 0
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
	
	** Generate Matrix of Regressors, shares, and outcome variable
	if ("`control_varlist'" ~= "") {
		mkmat `shiftshare_iv' `control_varlist' constant, matrix(Mn) //Matrix of IVs
		mkmat `endogenous_var' `control_varlist' constant, matrix(Gn)   //Matrix of regressors
	}
	else {
		mkmat `shiftshare_iv' constant, matrix(Mn)     //Matrix of IVs
		mkmat `endogenous_var' constant, matrix(Gn)    //Matrix of regressors
	}
	mkmat `share_varlist', matrix(ln)						//Matrix of Shares
	mkmat `dependant_var', matrix(tildeYn) 				//Dependent Variable   
	
	** OLS Estimates
	mat P1 = Mn' * Gn
	mat P2 = inv(Mn' * Mn)
	mat P3 = Mn' * tildeYn
	mat hat_theta = inv(P1'*P2*P1) * (P1'*P2*P3)
	mat e = tildeYn - Gn * hat_theta
	local hat_beta = hat_theta[1, 1]

	** Auxiliary variables
	mkmat `control_varlist' constant, matrix(tildeZn)
	mkmat `shiftshare_iv', matrix(tildeXn)
	mkmat `endogenous_var', matrix(tildeGn)
	local dim = rowsof(tildeXn)
	/*
	mat I_mat = I(`dim')
	mat Ydd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeYn
	mat Xdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeXn
	mat Gdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeGn
	mat Xddd = inv(ln'*ln) * (ln' * Xdd)
	*/
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
		mat R_raw = (e' * ln)'
		svmat R_raw, names(R_raw)
		svmat Xddd, names(Xddd)
		qui drop if mi(R_raw)
		
		gen R_raw_sq = R_raw^2
		gen Xddd_sq = Xddd^2

		mkmat R_raw_sq, matrix(R_sq)
		mkmat Xddd_sq, matrix(Xddd_sq)
		
		mat LambdaAKM = R_sq' * Xddd_sq

		mat variance = 1/`hatpi'^2 * inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
		local SE_AKM = sqrt(variance[1,1])
		local CI_low = `hat_beta' - `critical_value' * `SE_AKM'
		local CI_upp = `hat_beta' + `critical_value' * `SE_AKM'
		
		local tstat = `hat_beta' / `SE_AKM'
		local p = 2*(1 - normal(abs(`tstat')))
		
		** Output Results
		display " "
		display "The estimated coefficient is " %5.4f `hat_beta'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
		display "AKM               " %5.4f `SE_AKM'   "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp'
	}
	
	if `akmtype' == 0 {
		** Compute SE
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
		mat variance = 1/`hatpi'^2 * inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
		local SE_AKMnull_n = sqrt(variance[1,1])
		
		* Compute Confidence INterval
		local critical_value = invnormal(1-0.05/2)
		local critical2 = `critical_value'^2
		mat RY = Xdd' * Ydd
		mat RX = Xdd' * Gdd
		
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
		
		local SE_AKM = (`CI_upp' - `CI_low')/(2 * `critical_value')
		local CI_low = `hat_beta' - `critical_value' * `SE_AKM'
		local CI_upp = `hat_beta' + `critical_value' * `SE_AKM'
		
		local tstat = `hat_beta' / `SE_AKM'
		local p = 2*(1 - normal(abs(`tstat')))
		
		** Output Results
		display " "
		display "The estimated coefficient is " %5.4f `hat_beta'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
		display "AKM0              " %5.4f `SE_AKM'   "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp'
	}	
end
