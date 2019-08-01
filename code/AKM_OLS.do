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
	AKM0_nocluster
end

prog AKM_nocluster
	local control_var "t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource"

	*** (1) OLS regression
	** Full sample
	reg d_sh_empl d_tradeotch_pw_lag `control_var' [aw = weight]
	reg d_sh_empl d_tradeotch_pw_lag `control_var' [aw = weight], r
	reg d_sh_empl d_tradeotch_pw_lag `control_var' [aw = weight], cluster(state)

	* (AKM)
	set matsize 10000

	foreach var in `control_var' d_sh_empl d_tradeotch_pw_lag constant {
		replace `var' = `var' * sqrt(weight)
	}
	
	forvalues i=2(1)771 {
		replace var`i' = var`i' * sqrt(weight)
	}
	/*
	* Drop colinear share variables
	_rmcoll var*, force
	local share_emp_vars `r(varlist)'
	keep `r(varlist)'  d_sh_empl d_tradeotch_pw_lag `control_var' constant czone year weight
*/
	** Step (1) Regress Y_i on X_i and the controls Z_i
	reg d_sh_empl d_tradeotch_pw_lag `control_var' constant, noconstant
	local theta = _b[d_tradeotch_pw_lag]

	predict resid, residuals
	mkmat resid, matrix(resid_X)

	** Step (2) Regress X_i on Z_i
	qui reg d_tradeotch_pw_lag `control_var' constant, noconstant
	predict resid2, residuals
	mkmat resid2, matrix(resid_Xdot)

	** Step (3) Compute hat(X), the regression coefficients from regressing X_dot onto W
	qui reg resid2 var* 

	matrix X_hat = J(770, 1, 0)
	
	forvalues i = 2(1)771{
		loca j = `i'-1
		matrix X_hat[`j', 1] = _b[var`i']
	}
	/*
	local index = 1
	foreach var in `share_emp_vars' {
		capture confirm variable `var'
		
		if !_rc {
		   matrix X_hat[`index', 1] = _b[`var']
		   local index = `index' + 1
		}
	}
	*/

	mkmat var*, matrix(w)
	mat R = w' * resid_X
	mat R_sq = diag(R) * diag(R)

	mat X_dot_square = resid_Xdot' * resid_Xdot

	matrix C = J(770, 1, 0)
	forvalues i = 1(1)770 {
		mat C[`i', 1] = R[`i', 1]^2 * X_hat[`i', 1]^2
	}
	
	mat F = R_sq' * X_dot_square
	mat G = vecdiag(F)
	mata: st_matrix("H", rowsum(st_matrix("G")))
	local test = sqrt(H[1,1])
	display `test' 
	
	mata: st_matrix("A", colsum(st_matrix("C")))
	mata: st_matrix("B", rowsum(st_matrix("X_dot_square")))
	local se_beta = sqrt(A[1,1]) / B[1,1]

	local critical_value = invnormal(1-0.05/2)
	local ci_low = `theta' - `critical_value' * `se_beta'
	local ci_up  = `theta' + `critical_value' * `se_beta'

	display `theta'
	display `critical_value'
	display `se_beta' 
	display `ci_low' 
	display `ci_up'
	*/
end

prog AKM0_nocluster
	local control_var "t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource"

	*** (1) OLS regression
	** Full sample
	reg d_sh_empl d_tradeotch_pw_lag `control_var' [aw = weight]
	reg d_sh_empl d_tradeotch_pw_lag `control_var' [aw = weight], r
	reg d_sh_empl d_tradeotch_pw_lag `control_var' [aw = weight], cluster(state)

	* (AKM)
	set matsize 10000

	foreach var in `control_var' d_sh_empl d_tradeotch_pw_lag constant {
		replace `var' = `var' * sqrt(weight)
	}

	forvalues i=2(1)771 {
		replace var`i' = var`i' * sqrt(weight)
	}

	* Drop colinear share variables
	_rmcoll var*, force
	local share_emp_vars `r(varlist)'
	keep `r(varlist)'  d_sh_empl d_tradeotch_pw_lag `control_var' constant czone year weight
	

	** Step (1) Regress Y_i on X_i and the controls Z_i
	reg d_sh_empl d_tradeotch_pw_lag `control_var' constant, noconstant
	local theta = _b[d_tradeotch_pw_lag]
	predict resid, residual
	mkmat resid, matrix(resid_X)
	
	** Step (2) Regress Y_i - X_i beta_0 on Z_i
	local beta = 0
	gen restricted_dep_var = d_sh_empl - `beta' * d_tradeotch_pw_lag
	
	reg restricted_dep_var `control_var' constant, noconstant
	predict resid_restricted, residual
	mkmat resid_restricted, matrix(resid_restricted)

	** Step (3) Regress X_i on Z_i
	qui reg d_tradeotch_pw_lag `control_var' constant, noconstant
	predict resid2, residuals
	mkmat resid2, matrix(resid_Xdot)
	
	** Step (4) Compute hat(X), the regression coefficients from regressing X_dot onto W
	qui reg resid2 var* 

	matrix X_hat = J(770, 1, 0)
	local index = 1
	foreach var in `share_emp_vars' {
		capture confirm variable `var'
		
		if !_rc {
		   matrix X_hat[`index', 1] = _b[`var']
		   local index = `index' + 1
		}
	}

	mkmat var*, matrix(w)
	mat R = w' * resid_restricted
	mat R_s = w' * resid_X
	mat R_sq = diag(R) * diag(R)

	mat X_dot_square = resid_Xdot' * resid_Xdot

	matrix C = J(770, 1, 0)
	forvalues i = 1(1)770 {
		mat C[`i', 1] = R[`i', 1]^2 * X_hat[`i', 1]^2
	}

	mata: st_matrix("A", colsum(st_matrix("C")))
	mata: st_matrix("B", rowsum(st_matrix("X_dot_square")))
	local se_beta = sqrt(A[1,1]) / B[1,1]
	local B2 = B[1,1]
	display "se_beta = " `se_beta'
		
	local critical_value = invnormal(1-0.05/2)
	
	mat Q_minuspart_temp = w' * resid_Xdot
	mat Q_minuspart = J(770, 1, 0)
	forvalues i = 1(1)770 {
		mat Q_minuspart[`i', 1] = Q_minuspart_temp[`i', 1]^2 * X_hat[`i', 1]^2
	}
	
	mata: st_matrix("Q_temp", colsum(st_matrix("Q_minuspart")))
	local Q_t = Q_temp[1,1]

	
	local Q = `B2'^2  / `critical_value'^2 - `Q_t'
	display "Q = " `Q'
	
	mat A_temp = J(770, 1, 0)

	forvalues i = 1(1)770 {
		mat A_temp[`i', 1] = R_s[`i', 1] * X_hat[`i', 1]^2 * Q_minuspart_temp[`i', 1]
	} 
	mata: st_matrix("D", colsum(st_matrix("A_temp")))
	local A = D[1,1]/`Q'
	display `A'
	
	local ci_low = `theta' - `A' - sqrt(`A'^2 + `se_beta'^2/`Q' * `B2'^2)
	local ci_up  = `theta' - `A' + sqrt(`A'^2 + `se_beta'^2/`Q' * `B2'^2)
	display `theta'
	display `ci_low'
	display `ci_up'
end

main

cap log close
