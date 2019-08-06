program define AKM_OLS, eclass
	syntax, dependant_var(str) shiftshare_var(str) share_varlist(str) alpha(str) ///
			akmtype(str) [control_varlist(str) weight_var(str) beta0(str) path_cluster(str)]
	
	set more off
	set matsize 10000

	** OLS results without AKM adjustment
	local critical_value = invnormal(1-0.05/2)
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
	
	if "`beta0'" == "" {
		local beta0 = 0
	}
	
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
	mkmat `share_varlist', matrix(ln)						//Matrix of Shares
	mkmat `dependant_var', matrix(tildeYn) 				//Dependent Variable  
	
	** Estimate AKM standard error
	if `akmtype' == 1 {
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
			preserve
			
			** Get list of share variables by cluster
			use "`path_cluster'", clear
			qui tostring sector, replace

			qui levelsof sec_3d
			local sector_3d_list = r(levels)
			local num_cluster = 0

			foreach sector_3d in `sector_3d_list' {
				local num_cluster = `num_cluster' + 1
				qui levelsof sector if sec_3d == `sector_3d'
	
				local cluster_list`num_cluster' = r(levels)
			}
			
			restore
			
			mat LambdaAKM = J(1,1,0)
			forvalues i = 1(1)`num_cluster' {
				display "Compute SE of Sector " "`i'"
				local num_withincluster_category = 0
				
				foreach category in `cluster_list`i''{
					local num_withincluster_category = `num_withincluster_category' + 1
					mkmat emp_share`category' , matrix(temp`num_withincluster_category')
				} 
				
				matrix small_share_matrix = temp1
				
				* When having more than 1 categories within cluster
				if `num_withincluster_category' > 1 { 
					forvalues j = 2(1)`num_withincluster_category' {
						matrix small_share_matrix = small_share_matrix , temp`j'
					}
				}
				
				* Generate Xddd_cluster
				mat Xddd_cluster = J(`num_withincluster_category ', 1, 0)
				local temp_count = 0
				foreach category in `cluster_list`i''{
				local temp_count = `temp_count' + 1
					mat Xddd_cluster[`temp_count', 1] = Xddd[`category', 1]
				} 
				
				*
				mat RXcluster_raw = (e' * small_share_matrix)'
				svmat RXcluster_raw, names(RXcluster_raw)
				svmat Xddd_cluster, names(Xddd_cluster)
				
				preserve
				qui drop if mi(RXcluster)
				
				qui gen RXcluster = RXcluster_raw * Xddd_cluster
				mkmat RXcluster, matrix(RXcluster)
				
				mat RXcluster_sq = RXcluster * RXcluster'
				
				* Restore dataset and Drop variables
				restore
				drop RXcluster_raw* Xddd_cluster* RXcluster 
				
				mata : st_matrix("RXcluster_sq_rowsum", rowsum(st_matrix("RXcluster_sq")))
				mata : st_matrix("RXcluster_sq_fullsum", colsum(st_matrix("RXcluster_sq_rowsum")))
				mat LambdaAKM = LambdaAKM + RXcluster_sq_fullsum[1,1]
			}
		}
		
		mat variance = inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
		local SE_AKM = sqrt(variance[1,1])
		local CI_low = `coef' - `critical_value' * `SE_AKM'
		local CI_upp = `coef' + `critical_value' * `SE_AKM'
		
		local tstat = `coef' / `SE_AKM'
		local p = 2*(1 - normal(abs(`tstat')))
		
		** Output Results
		display " "
		display "The estimated coefficient is " %5.4f `coef'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
		display "AKM               " %5.4f `SE_AKM'   "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp'
	}
	
	** Estimate AKM0 standard error
	if `akmtype' == 0 {
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
		mat e_null = Ydd - Xdd * `beta0'
		
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
			preserve
			
			** Get list of share variables by cluster
			use "`path_cluster'", clear
			qui tostring sector, replace

			qui levelsof sec_3d
			local sector_3d_list = r(levels)
			local num_cluster = 0

			foreach sector_3d in `sector_3d_list' {
				local num_cluster = `num_cluster' + 1
				qui levelsof sector if sec_3d == `sector_3d'
	
				local cluster_list`num_cluster' = r(levels)
			}
			
			restore
			
			mat LambdaAKM = J(1,1,0)
			mat SXY = J(1, 1, 0)
			mat SYY = J(1, 1, 0)
			mat SXX = J(1, 1, 0)
			
			forvalues i = 1(1)`num_cluster' {
				display "Compute SE of Sector " "`i'"
				local num_withincluster_category = 0
				
				foreach category in `cluster_list`i''{
					local num_withincluster_category = `num_withincluster_category' + 1
					mkmat emp_share`category' , matrix(temp`num_withincluster_category')
				} 
				
				matrix small_share_matrix = temp1
				
				* When having more than 1 categories within cluster
				if `num_withincluster_category' > 1 { 
					forvalues j = 2(1)`num_withincluster_category' {
						matrix small_share_matrix = small_share_matrix , temp`j'
					}
				}
				
				* Generate Xddd_cluster
				mat Xddd_cluster = J(`num_withincluster_category ', 1, 0)
				local temp_count = 0
				foreach category in `cluster_list`i''{
				local temp_count = `temp_count' + 1
					mat Xddd_cluster[`temp_count', 1] = Xddd[`category', 1]
				} 
				
				*
				mat RXcluster_raw = (e_null' * small_share_matrix)'
				mat Ydd_lnCluster = (Ydd'* small_share_matrix)'
				mat Xdd_lnCluster = (Xdd'* small_share_matrix)'
				
				svmat RXcluster_raw, names(RXcluster_raw)
				svmat Xddd_cluster,  names(Xddd_cluster)
				svmat Ydd_lnCluster, names(Ydd_lnCluster)
				svmat Xdd_lnCluster, names(Xdd_lnCluster)
				
				preserve
				qui drop if mi(RXcluster)
				
				qui gen RXcluster  = RXcluster_raw * Xddd_cluster
				qui gen lnYCluster = Ydd_lnCluster * Xddd_cluster
				qui gen lnXCluster = Xdd_lnCluster * Xddd_cluster
				
				mkmat RXcluster, matrix(RXcluster)
				mkmat lnYCluster, matrix(lnYCluster)
				mkmat lnXCluster, matrix(lnXCluster)

				mat RXcluster_sq   = RXcluster * RXcluster'
				mat SXY_cluster_sq = lnXCluster * lnYCluster'
				mat SYY_cluster_sq = lnYCluster * lnYCluster'
				mat SXX_cluster_sq = lnXCluster * lnXCluster'
				
				* Restore dataset and Drop variables
				restore
				drop RXcluster_raw* Xddd_cluster* RXcluster Ydd_lnCluster* Xdd_lnCluster
				
				mata : st_matrix("RXcluster_sq_rowsum", rowsum(st_matrix("RXcluster_sq")))
				mata : st_matrix("RXcluster_sq_fullsum", colsum(st_matrix("RXcluster_sq_rowsum")))
				mata : st_matrix("SYY_cluster_sq_rowsum", rowsum(st_matrix("SYY_cluster_sq")))
				mata : st_matrix("SYY_cluster_sq_fullsum", colsum(st_matrix("SYY_cluster_sq_rowsum")))
				mata : st_matrix("SXX_cluster_sq_rowsum", rowsum(st_matrix("SXX_cluster_sq")))
				mata : st_matrix("SXX_cluster_sq_fullsum", colsum(st_matrix("SXX_cluster_sq_rowsum")))
				mata : st_matrix("SXY_cluster_sq_rowsum", rowsum(st_matrix("SXY_cluster_sq")))
				mata : st_matrix("SXY_cluster_sq_fullsum", colsum(st_matrix("SXY_cluster_sq_rowsum")))
				
				mat LambdaAKM = LambdaAKM + RXcluster_sq_fullsum[1,1]
				mat SXY = SXY + SXY_cluster_sq_fullsum[1,1]
				mat SYY = SYY + SYY_cluster_sq_fullsum[1,1]
				mat SXX = SXX + SXX_cluster_sq_fullsum[1,1]
			}
			
			preserve
			qui drop if mi(Xddd_sq) 
			mkmat Xddd_sq, matrix(Xddd_sq)
			restore
		}
		
		* Variance matrix
		mat variance = inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
		local SE_AKMnull_n = sqrt(variance[1,1])
		
		* Compute Confidence Interval
		local critical2 = `critical_value'^2
		mat RY = Xdd' * Ydd
		mat RX = Xdd' * Xdd
		
		if "`path_cluster'" == "" {
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
		
		** Output Results
		display " "
		display "The estimated coefficient is " %5.4f `coef'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r' 
		display "AKM0              " %5.4f `SE_AKM0'  "    " %5.4f `p' "    " %5.4f `CI_low' "    " %5.4f `CI_upp' 
		display "In AKM0 Estimation, the Confidence Interval Type is: Type " `CIType'
	}
end
