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
	
	* AKM_nocluster
	AKM0_nocluster
end

prog AKM_nocluster
	local control_var "t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource"
	
	set matsize 10000

	foreach var in `control_var' d_sh_empl d_tradeotch_pw_lag d_tradeusch_pw constant {
		replace `var' = `var' * sqrt(weight)
	}
	
	forvalues i = 1(1)792 {
		replace var`i' = var`i' * sqrt(weight)
	}
	
	* Drop colinear share variables
	_rmcoll var*, force
	local share_emp_vars `r(varlist)'
	keep `r(varlist)'  d_sh_empl d_tradeotch_pw_lag d_tradeusch_pw `control_var' constant czone year weight
	
	reg d_tradeusch_pw d_tradeotch_pw_lag `control_var' constant, noconstant
	ivregress 2sls d_sh_empl `control_var' constant (d_tradeusch_pw = d_tradeotch_pw_lag), noconstant
	
	** Generate Matrix of Regressors, shares, and outcome variable
	mkmat d_tradeotch_pw_lag `control_var' constant, matrix(Mn)   //Matrix of IVs
	mkmat d_tradeusch_pw `control_var' constant, matrix(Gn)   //Matrix of regressors
	mkmat var*, matrix(ln)								 //Matrix of Shares
	mkmat d_sh_empl, matrix(tildeYn) 					 //Dependent Variable  

	** OLS Estimates
	mat PM = Mn * inv(Mn' * Mn) * Mn'
	mat hat_theta = inv(Gn'*PM*Gn)*(Gn' * PM * tildeYn)
	mat e = tildeYn - Gn * hat_theta
	local hat_beta = hat_theta[1, 1]
    
	** Auxiliary variables
	mkmat `control_var' constant, matrix(tildeZn)
	mkmat d_tradeotch_pw_lag, matrix(tildeXn)
	mkmat d_tradeusch_pw, matrix(tildeGn)
	local dim = rowsof(tildeXn)
	
	mat I_mat = I(`dim')
	mat Ydd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeYn
	mat Xdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeXn
	mat Gdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeGn
	mat Xddd = inv(ln'*ln) * (ln' * Xdd)
	
	** First stage coefficient
	mat hat_thetaFS = inv(Mn'*Mn) * (Mn'*Gdd)
	local hatpi = hat_thetaFS[1,1]

	** Compute SE
	mat R_raw = (e' * ln)'
	svmat R_raw, names(R_raw)
	svmat Xddd, names(Xddd)
	drop if mi(R_raw)
	
	gen R_raw_sq = R_raw^2
	gen Xddd_sq = Xddd^2

	mkmat R_raw_sq, matrix(R_sq)
	mkmat Xddd_sq, matrix(Xddd_sq)
	
	mat LambdaAKM = R_sq' * Xddd_sq

	mat variance = 1/`hatpi'^2 * inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
	local SE = sqrt(variance[1,1])
	display "The Coefficient is: " %5.4f `hat_beta'
	display "The Standard Error is: " %5.4f `SE'
	
	local critical_value = invnormal(1-0.05/2)
	local CI_low = `hat_beta' - `critical_value' * `SE'
	local CI_upp = `hat_beta' + `critical_value' * `SE'
	display "The Confidence Interval is: [" %5.4f `CI_low' "," %5.4f `CI_upp' "]"
end


prog AKM0_nocluster
	local control_var "t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource"
	
	set matsize 10000

	foreach var in `control_var' d_sh_empl d_tradeotch_pw_lag d_tradeusch_pw constant {
		replace `var' = `var' * sqrt(weight)
	}
	
	forvalues i = 1(1)792 {
		replace var`i' = var`i' * sqrt(weight)
	}
	
	* Drop colinear share variables
	_rmcoll var*, force
	local share_emp_vars `r(varlist)'
	keep `r(varlist)'  d_sh_empl d_tradeotch_pw_lag d_tradeusch_pw `control_var' constant czone year weight
	
	reg d_tradeusch_pw d_tradeotch_pw_lag `control_var' constant, noconstant
	ivregress 2sls d_sh_empl `control_var' constant (d_tradeusch_pw = d_tradeotch_pw_lag), noconstant
	
	** Generate Matrix of Regressors, shares, and outcome variable
	mkmat d_tradeotch_pw_lag `control_var' constant, matrix(Mn)   //Matrix of IVs
	mkmat d_tradeusch_pw `control_var' constant, matrix(Gn)   //Matrix of regressors
	mkmat var*, matrix(ln)								 //Matrix of Shares
	mkmat d_sh_empl, matrix(tildeYn) 					 //Dependent Variable  

	** OLS Estimates
	mat PM = Mn * inv(Mn' * Mn) * Mn'
	mat hat_theta = inv(Gn'*PM*Gn)*(Gn' * PM * tildeYn)
	mat e = tildeYn - Gn * hat_theta
	local hat_beta = hat_theta[1, 1]
    
	** Auxiliary variables
	mkmat `control_var' constant, matrix(tildeZn)
	mkmat d_tradeotch_pw_lag, matrix(tildeXn)
	mkmat d_tradeusch_pw, matrix(tildeGn)
	local dim = rowsof(tildeXn)
	
	mat I_mat = I(`dim')
	mat Ydd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeYn
	mat Xdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeXn
	mat Gdd = (I_mat - tildeZn * inv(tildeZn'*tildeZn) * tildeZn') * tildeGn
	mat Xddd = inv(ln'*ln) * (ln' * Xdd)
	
	** First stage coefficient
	mat hat_thetaFS = inv(Mn'*Mn) * (Mn'*Gdd)
	local hatpi = hat_thetaFS[1,1]
	
	** Compute SE
	local beta0 = 0
	mat e_null = Ydd - Xdd * `beta0'
	
	mat R_raw = (e_null' * ln)'
	svmat R_raw, names(R_raw)
	svmat Xddd, names(Xddd)
	drop if mi(R_raw)
	
	gen R_raw_sq = R_raw^2
	gen Xddd_sq = Xddd^2

	mkmat R_raw_sq, matrix(R_sq)
	mkmat Xddd_sq, matrix(Xddd_sq)
	
	mat LambdaAKM = R_sq' * Xddd_sq

    * Variance matrix
	mat variance = 1/`hatpi'^2 * inv(Xdd'*Xdd) * LambdaAKM * inv(Xdd'*Xdd)
	local SE_AKMnull_n = sqrt(variance[1,1])
	display `SE_AKMnull_n'
	
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
	display `Q'
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
    
	local SE = (`CI_upp' - `CI_low')/(2 * `critical_value')

	display "The Standard Error is: " %5.4f `SE'
	display "The Confidence Interval is: [" %5.4f `CI_low' "," %5.4f `CI_upp' "]"
end

main

cap log close
