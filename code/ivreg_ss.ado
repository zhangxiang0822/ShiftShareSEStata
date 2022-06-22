program define ivreg_ss, eclass
	version 14.0
	syntax varlist(numeric min=1 max=1), endogenous_var(varlist numeric min=1 max=1) ///
		   shiftshare_iv(varlist numeric min=1 max=1) share_varlist(varlist numeric min=1) ///
		   [control_varlist(varlist numeric min=1) weight_var(varlist numeric min=1 max=1) ///
		    akmtype(str) beta0(real 0.0) alpha(str) path_cluster(str) cluster_var(str) firststage(integer 0)]
	
	set more off
	set matsize 10000
	
	preserve
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
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
						 alpha(`alpha') control_varlist(`control_varlist') weight_var(`weight_var') akmtype(`akmtype') ///
						 path_cluster(`path_cluster') cluster_var(`cluster_var')
				}
				else {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
						 alpha(`alpha') control_varlist(`control_varlist') weight_var(`weight_var') akmtype(`akmtype')
				}
			}
			else {
				if "`cluster_var'" ~= "" {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') control_varlist(`control_varlist') akmtype(`akmtype') ///
							 path_cluster(`path_cluster') cluster_var(`cluster_var') 
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
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') weight_var(`weight_var') akmtype(`akmtype') ///
							 path_cluster(`path_cluster') cluster_var(`cluster_var')
				} 
				else {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') weight_var(`weight_var') akmtype(`akmtype')
				}
			}
			else {
				if "`cluster_var'" ~= "" {
					reg_ss `endogenous_var', shiftshare_var(`shiftshare_iv') share_varlist(`share_varlist') ///
							 alpha(`alpha') akmtype(`akmtype') path_cluster(`path_cluster') cluster_var(`cluster_var')
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
		local p_homo  = 2 * (1 - normal(abs(`z_homo')))
		local CI_low_homo = _b[`endogenous_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`endogenous_var'] + `critical_value' * `SE_homo'
		
		qui ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], r 
		local SE_r = _se[`endogenous_var']
		local z_r  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_r  = 2 * (1 - normal(abs(`z_r')))
		local CI_low_r = _b[`endogenous_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`endogenous_var'] + `critical_value' * `SE_r'
		
		* ivregress 2sls `dependant_var' `control_varlist' (`endogenous_var' = `shiftshare_iv') [aw = `weight_var'], cluster(state)
	}
	else {
		qui ivregress 2sls `dependant_var' `control_varlist'  (`endogenous_var' = `shiftshare_iv') 
		local SE_homo = _se[`endogenous_var']
		local z_homo  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_homo  = 2 * (1 - normal(abs(`z_homo')))
		local CI_low_homo = _b[`endogenous_var'] - `critical_value' * `SE_homo'
		local CI_upp_homo = _b[`endogenous_var'] + `critical_value' * `SE_homo'
		
		qui ivregress 2sls `dependant_var' `control_varlist'  (`endogenous_var' = `shiftshare_iv'), vce(robust)
		local SE_r = _se[`endogenous_var']
		
		local z_r  = _b[`endogenous_var']/_se[`endogenous_var']
		local p_r  = 2 * (1 - normal(abs(`z_r')))
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
		display "Warning: You have collinear share variables (Share matrix has colinear columns)"
	}
	
	** Generate Matrix of Regressors, shares, and outcome variable
	if ("`control_varlist'" ~= "") {
		local IV_varlist `shiftshare_iv' `controls'
		local regressor_varlist `endogenous_var' `controls'
	}
	else {
		local IV_varlist `shiftshare_iv' constant
		local regressor_varlist `endogenous_var' constant
	}
	
	marksample touse
	markout `touse' `shiftshare_var' `regressor_varlist' `share_varlist', strok
	
	** Estimate AKM standard error
	if `akmtype' == 1 {
	
		tempname ln Xdd Xddd e 
		mata: s_AKM1_consvec_iv("`IV_varlist'", "`regressor_varlist'", "`share_varlist'", "`dependant_var'", "`touse'")
		
		** commonly used matrices
		mat `ln'	= r(ln)
		mat `Xdd'	= r(Xdd)
		mat `Xddd'	= r(Xddd)
		mat `e'		= r(e)
		
		** Get regression coefficient
		local coef = `r(b)'
		local fs_coef = `r(fs_b)'
		
		** Check number of regions and number of sectors
		local num_counties = rowsof(`ln')
		local num_sector1  = colsof(`ln')
	
		if (`num_counties' < `num_sector1') {
			display "ERROR: You have more number of sectors than regions"
			exit
		}
		
		restore
		
		** Compute SE
		if "`path_cluster'" == "" {
			mata: s_AKM1_nocluster_iv("`e'", "`ln'", "`Xdd'", "`Xddd'", `coef', `fs_coef', `critical_value')
		}
		else {		
			preserve
			use "`path_cluster'", clear
				
			local 0 `cluster_var'
			syntax varlist
			
			marksample alluse
			markout `alluse' `cluster_var', strok
			
			qui sum `cluster_var'
			local num_sector_obs = r(N)
			
			** translate sector data into matrix
			tempname sec_vec_full sec_vec_unique
			
			mata: s_readsec_iv("`cluster_var'", "`alluse'")
			mat `sec_vec_full' = r(sec_vec)
			
			** get unique values of `cluster_var'
			qui duplicates drop `cluster_var', force
			mata: s_readsec_iv("`cluster_var'", "`alluse'")
			mat `sec_vec_unique' = r(sec_vec)
			
			** Compute SE
			mata: s_AKM1_cluster_iv("`e'", "`ln'", "`Xdd'", "`Xddd'", `beta0', `critical_value', "`sec_vec_full'", "`sec_vec_unique'", `coef', `fs_coef')
			
			restore
		}
		
		** Display results
		display " "
		display "The estimated coefficient is " %5.4f `coef'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r'
		display "AKM               " %5.4f `r(se)'  "    " %5.4f `r(p)' "    " %5.4f `r(CIl)' "    " %5.4f `r(CIu)'
	}
	
	** Estimate AKM0 standard error
	if `akmtype' == 0 {
	
		tempname ln Xdd Xddd Ydd Gdd e e_null 
		mata: s_AKM0_consvec_iv("`IV_varlist'", "`regressor_varlist'", "`share_varlist'", "`dependant_var'", "`touse'", `beta0')
			
		** commonly used matrices
		mat `ln'		= r(ln)
		mat `Xdd'		= r(Xdd)
		mat `Xddd'		= r(Xddd)
		mat `Ydd'		= r(Ydd)
		mat `Gdd'		= r(Gdd)
		mat `e'			= r(e)
		mat `e_null'	= r(e_null)
		
		** Get regression coefficient
		local coef = `r(b)'
		local fs_coef = `r(fs_b)'
		
		** Check number of regions and number of sectors
		local num_counties = rowsof(`ln')
		local num_sector1  = colsof(`ln')
	
		if (`num_counties' < `num_sector1') {
			display "ERROR: You have more number of sectors than regions"
			exit
		}
		
		restore
		
		** Compute standard error
		if "`path_cluster'" == "" {
			mata: s_AKM0_nocluster_iv("`e_null'", "`ln'", "`Xdd'", "`Xddd'", "`Ydd'", "`Gdd'", `coef', `fs_coef', `beta0', `critical_value')
		}
		else {
			preserve
			use "`path_cluster'", clear
				
			local 0 `cluster_var'
			syntax varlist
			
			marksample alluse
			markout `alluse' `cluster_var', strok
			
			qui sum `cluster_var'
			local num_sector_obs = r(N)
			
			** translate sector data into matrix
			tempname sec_vec_full sec_vec_unique
			
			mata: s_readsec_iv("`cluster_var'", "`alluse'")
			mat `sec_vec_full' = r(sec_vec)
			
			** get unique values of `cluster_var'
			qui duplicates drop `cluster_var', force
			mata: s_readsec_iv("`cluster_var'", "`alluse'")
			mat `sec_vec_unique' = r(sec_vec)
			
			** Compute SE
			mata: s_AKM0_cluster_iv("`e_null'", "`ln'", "`Xdd'", "`Xddd'", "`Ydd'", "`Gdd'", `coef', `fs_coef', `beta0', `critical_value', "`sec_vec_full'", "`sec_vec_unique'")
			
			restore
		}
		
		** Display results
		display " "
		display "The estimated coefficient is " %5.4f `coef'
		display "Inference"
		display "               Std. Error   p-value   Lower CI   Upper CI"
		display "Homoscedastic     " %5.4f `SE_homo'  "    " %5.4f `p_homo' "    " %5.4f `CI_low_homo' "    " %5.4f `CI_upp_homo' 
		display "EHW               " %5.4f `SE_r'     "    " %5.4f `p_r' "    " %5.4f `CI_low_r' "    " %5.4f `CI_upp_r'
		
		if `r(CItype)' == 1 {
			display "AKM0              " %5.4f `r(se)'  "    " %5.4f `r(p)' "    " %5.4f `r(CIl)' "    " %5.4f `r(CIu)'
		}
		else if `r(CItype)' == 2 {
			display "Inference"
			display "               Std. Error   p-value   CI"
			display "AKM0              " %5.4f `r(se)'  "    " %5.4f `r(p)' "    " "=(-Inf, " %5.4f `r(CIl)' "] + [" %5.4f `r(CIu)' " , Inf)"
		}
		else {
			display "Inference"
			display "               Std. Error   p-value   CI"
			display "AKM0              " %5.4f `r(se)'  "    " %5.4f `r(p)' "    " "=(-Inf, Inf)"
		}
	}
	
	** Save results for return
	ereturn scalar b = `coef'
	ereturn local se = `r(se)'
	ereturn local CI_upp = `r(CIu)'
	ereturn local CI_low = `r(CIl)'
	ereturn scalar p = `r(p)'
	ereturn scalar tstat = `r(tstat)'
	
	if `firststage' != 0 {
		ereturn local se_firststage = `se_firststage'
		ereturn local b_firststage  = `b_firststage'
		ereturn local CI_upp_firststage = `CI_upp_firststage'
		ereturn local CI_low_firststage = `CI_low_firststage'
		ereturn local p_firststage = `p_firststage'
		ereturn scalar tstat = `tstat_firststage'
	}
	
end

*******************************************************************************
*************************** BEGIN MATA CODE ***********************************
*******************************************************************************

version 14.0
set matastrict on

mata:

void s_AKM1_consvec_iv(string scalar IV_var,
					   string scalar regressor_var, 
					   string scalar share_varlist, 
					   string scalar dependant_var,
					   string scalar touse)
{
	real matrix Mn
	real matrix Gn
	real matrix ln
	real colvector tildeYn
	real matrix P1
	real matrix P2
	real matrix P3
	real matrix hat_theta
	real matrix hat_beta
	real scalar coef
	real matrix tildeZn
	real matrix tildeXn
	real matrix A
	real matrix Xdd
	real matrix Xddd
	real matrix Ydd
	real matrix Gdd
	real matrix e_ln
	
	// Generate Matrix of Regressors, shares, and outcome variable
	//Mn: Matrix of IVs
	//Gn: Matrix of regressors
	//ln: Matrix of Shares
	//tildeYn: Vector of Dependent Variable
	
	st_view(Mn=., ., tokens(IV_var), touse) 	
	st_view(Gn=., ., tokens(regressor_var), touse) 	
	st_view(ln=., ., tokens(share_varlist), touse) 		
	st_view(tildeYn=., ., tokens(dependant_var), touse) 
	
	P1 = Mn' * Gn
	P2 = invsym(Mn' * Mn)
	P3 = Mn' * tildeYn
	hat_theta = invsym(P1' * P2 * P1) * (P1' * P2 * P3)

	e = tildeYn - Gn*hat_theta
	hat_beta = hat_theta[1, 1]

	tildeZn = Mn[|1,2\.,.|]
	tildeXn = Mn[|1,1\.,1|]
	
	A = tildeZn * invsym(tildeZn' * tildeZn)
	Xdd = tildeXn - A * (tildeZn' * tildeXn)
	Ydd = tildeYn - A * (tildeZn' * tildeYn)
	Gdd = Gn[|1,1\.,1|] - A * (tildeZn' * Gn[|1,1\.,1|])
	Xddd = invsym(ln' *ln) * (ln' * Xdd)
	
	// First stage coefficient
	hat_thetaFS = invsym(Mn' * Mn) * (Mn' * Gdd)
	hatpi = hat_thetaFS[1,1]
	
	// convert from view to matrix that can be returned
	ln_q   = ln[|.,.\.,.|]
	Xdd_q  = Xdd[|.,.\.,.|]
	Xddd_q = Xddd[|.,.\.,.|]
	e_q    = e[|.,.\.,.|]
	
	st_matrix("r(ln)", ln_q)
	st_matrix("r(Xdd)", Xdd_q)
	st_matrix("r(Xddd)", Xddd_q)
	st_matrix("r(e)", e_q)
	st_numscalar("r(b)", hat_beta)
	st_numscalar("r(fs_b)", hatpi)
}	

void s_AKM0_consvec_iv(string scalar IV_var,
					   string scalar regressor_var, 
					   string scalar share_varlist, 
					   string scalar dependant_var,
					   string scalar touse,
					   scalar beta0)
{
	real matrix Mn
	real matrix Gn
	real matrix ln
	real colvector tildeYn
	real matrix P1
	real matrix P2
	real matrix P3
	real matrix hat_theta
	real matrix hat_beta
	real scalar coef
	real matrix tildeZn
	real matrix tildeXn
	real matrix A
	real matrix Xdd
	real matrix Xddd
	real matrix Ydd
	real matrix Gdd
	real matrix e_ln
	real matrix e_null
	
	//Generate Matrix of Regressors, shares, and outcome variable
	//Mn: Matrix of IVs
	//Gn: Matrix of regressors
	//ln: Matrix of Shares
	//tildeYn: Vector of Dependent Variable
	
	st_view(Mn=., ., tokens(IV_var), touse) 	
	st_view(Gn=., ., tokens(regressor_var), touse) 	
	st_view(ln=., ., tokens(share_varlist), touse) 		
	st_view(tildeYn=., ., tokens(dependant_var), touse) 
	
	P1 = Mn' * Gn
	P2 = invsym(Mn' * Mn)
	P3 = Mn' * tildeYn
	hat_theta = invsym(P1' * P2 * P1) * (P1' * P2 * P3)

	e = tildeYn - Gn*hat_theta
	
	hat_beta = hat_theta[1, 1]

	tildeZn = Mn[|1,2\.,.|]
	tildeXn = Mn[|1,1\.,1|]
	
	A = tildeZn * invsym(tildeZn' * tildeZn)
	Xdd = tildeXn - A * (tildeZn' * tildeXn)
	Ydd = tildeYn - A * (tildeZn' * tildeYn)
	Gdd = Gn[|1,1\.,1|] - A * (tildeZn' * Gn[|1,1\.,1|])
	Xddd = invsym(ln' *ln) * (ln' * Xdd)
	e_null = Ydd - Gdd * beta0
	
	// First stage coefficient
	hat_thetaFS = invsym(Mn' * Mn) * (Mn' * Gdd)
	hatpi = hat_thetaFS[1,1]
	
	// convert from view to matrix that can be returned
	ln_q   = ln[|.,.\.,.|]
	Xdd_q  = Xdd[|.,.\.,.|]
	Xddd_q = Xddd[|.,.\.,.|]
	e_q    = e[|.,.\.,.|]
	Gdd_q  = Gdd[|.,.\.,.|]
	Ydd_q  = Ydd[|.,.\.,.|]
	e_null_q = e_null[|.,.\.,.|]
	
	st_matrix("r(ln)", ln_q)
	st_matrix("r(Xdd)", Xdd_q)
	st_matrix("r(Xddd)", Xddd_q)
	st_matrix("r(e)", e_q)
	st_matrix("r(Ydd)", Ydd_q)
	st_matrix("r(Gdd)", Gdd_q)
	st_matrix("r(e_null)", e_null_q)
	st_numscalar("r(b)", hat_beta)
	st_numscalar("r(fs_b)", hatpi)
}	



void s_readsec_iv(string scalar cluster_var, string scalar touse) 
{
	real matrix sec_vec
	
	st_view(sec_vec=., ., tokens(cluster_var), touse)
	
	sec_vec_q = sec_vec[|.,.\.,.|]
	
	st_matrix("r(sec_vec)", sec_vec_q)
}

void s_AKM1_nocluster_iv(string scalar e_matrix,	 
					  string scalar ln_matrix,   
					  string scalar Xdd_matrix,  
					  string scalar Xddd_matrix, 
					  scalar coef, 	
					  scalar fs_coef,
					  scalar critical_value)
{	
	real matrix e
	real matrix ln
	real matrix Xdd
	real matrix Xddd

	// Read in data
	e    = st_matrix(e_matrix)
	ln   = st_matrix(ln_matrix)
	Xdd  = st_matrix(Xdd_matrix)
	Xddd = st_matrix(Xddd_matrix)
	
	// Compute variance matrix
	R = (e' * ln) :* (e' * ln)
	LambdaAKM = R * (Xddd :* Xddd)
	
	variance = (1/fs_coef^2) * invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
		
	//AKM se, p-value, and t-stat
	se_AKM = sqrt(variance[1,1])
	tstat = coef / se_AKM
	p_value = 2 * (1 - normal(abs(tstat)))
	
	CIl = coef - critical_value * se_AKM
	CIu = coef + critical_value * se_AKM
	
	st_numscalar("r(se)", se_AKM)
	st_numscalar("r(tstat)", tstat)
	st_numscalar("r(p)", p_value)
	st_numscalar("r(CIl)", CIl)
	st_numscalar("r(CIu)", CIu)
}

void s_AKM0_nocluster_iv(string scalar e_null_matrix,	 
					  string scalar ln_matrix,   
					  string scalar Xdd_matrix,  
					  string scalar Xddd_matrix, 
					  string scalar Ydd_matrix,
					  string scalar Gdd_matrix,
					  scalar coef, 		
					  scalar fs_coef,
					  scalar beta0,
					  scalar critical_value)
{	
	real matrix e_null
	real matrix ln
	real matrix Xdd
	real matrix Xddd
	real matrix Ydd
	real matrix Gdd
	real scalar critical2
	real scalar RY
	real scalar RX
	real matrix lnY
	real matrix lnX
	real scalar SXY
	real scalar SXX
	real scalar SYY
	real scalar Q
	real scalar Delta
	real scalar CIu
	real scalar CIl
	
	// Read in data
	e_null  = st_matrix(e_null_matrix)
	ln   	= st_matrix(ln_matrix)
	Xdd  	= st_matrix(Xdd_matrix)
	Xddd 	= st_matrix(Xddd_matrix)
	Ydd  	= st_matrix(Ydd_matrix)
	Gdd  	= st_matrix(Gdd_matrix)
	
	// Compute variance matrix
	R = (e_null' * ln) :* (e_null' * ln)
	LambdaAKM = R * (Xddd :* Xddd)
	
	variance = (1/fs_coef^2) * invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
	
	se_AKMnull = sqrt(variance[1,1])
	tstat = (coef - beta0) / se_AKMnull
	
	p_value = 2 * (1 - normal(abs(tstat)))
	
	// Compute confidence interval
	critical2 = critical_value^2
    RY = Xdd' * Ydd
    RX = Xdd' * Gdd
	
	lnY =  Ydd' * ln 
    lnX =  Gdd' * ln
	
	XX = Xddd :* Xddd
    SXY = (lnY :* lnX) * XX
    SXX = (lnX :* lnX) * XX
    SYY = (lnY :* lnY) * XX
	
	Q = (RX^2)  /critical2 - SXX
    Delta = (RY * RX - critical2 * SXY)^2 - (RX^2 - critical2 * SXX)*(RY^2 - critical2 * SYY)
	
	if (Q > 0) {
			CIl = ( (RY*RX - critical2*SXY) - Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CIu = ( (RY*RX - critical2*SXY) + Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CItype = 1
		}
	else {
		if (Delta > 0) {
			CIl = ( (RY*RX - critical2*SXY) + Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CIu = ( (RY*RX - critical2*SXY) - Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CItype = 2
		}
		else {
			CIl = -10000000
			CIu = 10000000
			CItype = 3
		}
	}
	
	SE_AKM0 = (CIu - CIl)/(2 * critical_value)
		
	st_numscalar("r(se)", SE_AKM0)
	st_numscalar("r(tstat)", tstat)
	st_numscalar("r(p)", p_value)
	st_numscalar("r(CIl)", CIl)
	st_numscalar("r(CIu)", CIu)
	st_numscalar("r(CItype)", CItype)
}

void s_AKM1_cluster_iv(string scalar e_matrix,	 			
					string scalar ln_matrix,   			
					string scalar Xdd_matrix,  			
					string scalar Xddd_matrix, 			
					scalar beta0, 						
					scalar critical_value, 				
					string scalar sec_vec_full_matrix, 	
					string scalar sec_vec_unique_matrix, 
					scalar coef,
					scalar fs_coef
					)
{	
	real matrix e
	real matrix ln
	real matrix Xdd
	real matrix Xddd
	real matrix sec_vec_full
	real matrix sec_vec_unique
	real scalar LambdaAKM
	real matrix RXcluster
	real matrix lncluster
	
	// Read in data
	e    = st_matrix(e_matrix)
	ln   = st_matrix(ln_matrix)
	Xdd  = st_matrix(Xdd_matrix)
	Xddd = st_matrix(Xddd_matrix)
	sec_vec_full = st_matrix(sec_vec_full_matrix)
	sec_vec_unique = st_matrix(sec_vec_unique_matrix)
	
	// Compute variance matrix
	LambdaAKM  = 0
	
	startcol = 1
	endcol = 1
	
	nrow_full = rows(sec_vec_full)
	nrow_unique = rows(sec_vec_unique)
	norder = 1
	
	for (i = 1; i <= nrow_unique; i++) {
		sector = sec_vec_unique[i, 1]
		
		flag = 0
		
		for (j = 1; j <= nrow_full; j++) {
			if (sec_vec_full[j, 1] == sector) {
				if (flag == 0) {
					select_vec = (j)
					flag = 1
				}
				else {
					select_vec = (select_vec \ j)
				}
			}
		}
		
		lnCluster = ln[. , select_vec]
		XdddCluster = Xddd[select_vec, .]		
		RXCluster = (( e'*lnCluster  )') :* XdddCluster
		
		LambdaAKM  = LambdaAKM  + sum(RXCluster * RXCluster')
	}
	
	variance = (1/fs_coef^2) * invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
	
	//AKM se, p-value, and t-stat
	se_AKM = sqrt(variance[1,1])
	tstat = (coef - beta0) / se_AKM
	
	p_value = 2 * (1 - normal(abs(tstat)))
	
	CIl = coef - critical_value * se_AKM
	CIu = coef + critical_value * se_AKM
	
	st_numscalar("r(se)", se_AKM)
	st_numscalar("r(tstat)", tstat)
	st_numscalar("r(p)", p_value)
	st_numscalar("r(CIl)", CIl)
	st_numscalar("r(CIu)", CIu)
}

void s_AKM0_cluster_iv(string scalar e_null_matrix,	 
					   string scalar ln_matrix,   
					   string scalar Xdd_matrix,  
					   string scalar Xddd_matrix, 
					   string scalar Ydd_matrix,
					   string scalar Gdd_matrix,
					   scalar coef, 		
					   scalar fs_coef,
					   scalar beta0,
					   scalar critical_value,
					   string scalar sec_vec_full_matrix, 	
					   string scalar sec_vec_unique_matrix) 
{	
	real matrix e_null
	real matrix ln
	real matrix Xdd
	real matrix Xddd
	real matrix Ydd
	real matrix sec_vec_full
	real matrix sec_vec_unique
	real scalar LambdaAKM
	real scalar SYY
	real scalar SXY
	real scalar SXX
	real matrix RXCluster
	real matrix lnCluster
	real matrix lnYCluster
	real matrix lnXCluster
	real scalar critical2
	real scalar RY
	real scalar RX
	real matrix lnY
	real matrix lnX
	real scalar Q
	real scalar Delta
	real scalar CIu
	real scalar CIl
	
	// Read in data
	e_null = st_matrix(e_null_matrix)
	ln     = st_matrix(ln_matrix)
	Xdd    = st_matrix(Xdd_matrix)
	Xddd   = st_matrix(Xddd_matrix)
	Ydd    = st_matrix(Ydd_matrix)
	Gdd    = st_matrix(Gdd_matrix)
	sec_vec_full = st_matrix(sec_vec_full_matrix)
	sec_vec_unique = st_matrix(sec_vec_unique_matrix)
	
	// Compute variance matrix
	LambdaAKM  = 0
	SYY = 0
	SXY = 0
	SXX = 0
	
	startcol = 1
	endcol = 1
	
	nrow_full = rows(sec_vec_full)
	nrow_unique = rows(sec_vec_unique)
	norder = 1
	
	for (i = 1; i <= nrow_unique; i++) {
		sector = sec_vec_unique[i, 1]
		
		flag = 0
		
		for (j = 1; j <= nrow_full; j++) {
			if (sec_vec_full[j, 1] == sector) {
				if (flag == 0) {
					select_vec = (j)
					flag = 1
				}
				else {
					select_vec = (select_vec \ j)
				}
			}
		}
		
		lnCluster = ln[. , select_vec]
		XdddCluster = Xddd[select_vec, .]
		
		RXCluster = (( e_null'*lnCluster  )') :* XdddCluster
		lnYCluster = (( Ydd'*lnCluster  )') :* XdddCluster
        lnXCluster = (( Gdd'*lnCluster  )') :* XdddCluster
		
		LambdaAKM  = LambdaAKM  + sum(RXCluster * RXCluster')
		SXY = SXY + sum(lnYCluster * lnXCluster')
        SXX = SXX + sum(lnXCluster * lnXCluster')
        SYY = SYY + sum(lnYCluster * lnYCluster')
	}
	
	variance = (1/fs_coef^2) * invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
	
	// AKM se, p-value, and t-stat
	se_AKM = sqrt(variance[1,1])
	tstat = (coef - beta0) / se_AKM
	
	p_value = 2 * (1 - normal(abs(tstat)))
	
	// Compute confidence interval
	critical2 = critical_value^2
    RY = Xdd' * Ydd
    RX = Xdd' * Gdd
	
	Q = (RX^2)  /critical2 - SXX
    Delta = (RY * RX - critical2 * SXY)^2 - (RX^2 - critical2 * SXX)*(RY^2 - critical2 * SYY)
	
	if (Q > 0) {
			CIl = ( (RY*RX - critical2*SXY) - Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CIu = ( (RY*RX - critical2*SXY) + Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CItype = 1
		}
	else {
		if (Delta > 0) {
			CIl = ( (RY*RX - critical2*SXY) + Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CIu = ( (RY*RX - critical2*SXY) - Delta^(1/2) )/(RX^2 - critical2 * SXX)
			CItype = 2
		}
		else {
			CIl = -10000000
			CIu = 10000000
			CItype = 3
		}
	}
	
	SE_AKM0 = (CIu - CIl)/(2 * critical_value)
		
	st_numscalar("r(se)", SE_AKM0)
	st_numscalar("r(tstat)", tstat)
	st_numscalar("r(p)", p_value)
	st_numscalar("r(CIl)", CIl)
	st_numscalar("r(CIu)", CIu)
	st_numscalar("r(CItype)", CItype)
}

end		// end mata section
