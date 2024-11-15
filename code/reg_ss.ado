program define reg_ss, eclass
	version 14.0
	syntax varlist(numeric min=1 max=1), ///
		   shiftshare_var(varlist numeric min=1 max=1) share_varlist(varlist numeric min=1) ///
		   [control_varlist(varlist numeric min=1) weight_var(varlist numeric min=1 max=1) ///
		    akmtype(str) beta0(real 0.0) alpha(str) path_cluster(str) cluster_var(str)]
			
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
	
	** Some locals
	local critical_value = invnormal(1-`alpha'/2)
	
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
		local p_r  = 2*ttail(e(df_r),abs(`t_r'))
		local CI_low_r = _b[`shiftshare_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`shiftshare_var'] + `critical_value' * `SE_r'
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
		local p_r  = 2*ttail(e(df_r),abs(`t_r'))
		local CI_low_r = _b[`shiftshare_var'] - `critical_value' * `SE_r'
		local CI_upp_r = _b[`shiftshare_var'] + `critical_value' * `SE_r'
		
		* reg `dependant_var' `shiftshare_var' `control_varlist', cluster(state)
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
	
	capture _rmcoll `control_varlist' constant, force nocons
	if _rc == 0{
		_rmcoll `control_varlist' constant, force nocons
		local controls `r(varlist)'
	}
	else {
		display "You must manually generate dummy variables, instead of using i.XXX"
		exit
	}
	
	_rmcoll `share_varlist', force
	local num_omit_var = `r(k_omitted)'
	if `num_omit_var' > 0 {
		display "Warning: The share matrix has collinear columns. Collinear sectoral shocks have been dropped"
	}
	 
	** Generate Matrix of Regressors, shares, and outcome variable
	local regressor_var `shiftshare_var' `controls'
	
	marksample touse
	markout `touse' `shiftshare_var' `regressor_var' `share_varlist', strok
	
	** Estimate AKM standard error
	if `akmtype' == 1 {
	
		tempname ln Xdd Xddd e e_ln
		mata: s_AKM1_consvec("`regressor_var'", "`share_varlist'", "`dependant_var'", "`touse'")
			
		** commonly used matrices
		mat `ln'		= r(ln)
		mat `Xdd'		= r(Xdd)
		mat `Xddd'		= r(Xddd)
		mat `e'			= r(e)
		mat `e_ln'		= r(e_ln)
		
		** Get regression coefficient
		local coef = `r(b)'
		
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
		
			mata: s_AKM1_nocluster("`e'", "`ln'", "`Xdd'", "`Xddd'", `r(b)', `critical_value')
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
			
			mata: s_readsec("`cluster_var'", "`alluse'")
			mat `sec_vec_full' = r(sec_vec)
			
			** get unique values of `cluster_var'
			qui duplicates drop `cluster_var', force
			mata: s_readsec("`cluster_var'", "`alluse'")
			mat `sec_vec_unique' = r(sec_vec)
			
			** Compute SE
			mata: s_AKM1_cluster("`e'", "`ln'", "`Xdd'", "`Xddd'", `beta0', `critical_value', "`sec_vec_full'", "`sec_vec_unique'", `coef')
			
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
	
		tempname ln Xdd Xddd Ydd e e_ln e_null
		mata: s_AKM0_consvec("`regressor_var'", "`share_varlist'", "`dependant_var'", "`touse'", `beta0')
			
		** commonly used matrices
		mat `ln'		= r(ln)
		mat `Xdd'		= r(Xdd)
		mat `Xddd'		= r(Xddd)
		mat `Ydd'		= r(Ydd)
		mat `e'			= r(e)
		mat `e_ln'		= r(e_ln)
		mat `e_null'	= r(e_null)
		
		** Get regression coefficient
		local coef = `r(b)'
		
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
			
			mata: s_AKM0_nocluster("`e_null'", "`ln'", "`Xdd'", "`Xddd'", "`Ydd'", `coef', `beta0', `critical_value')
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
			
			mata: s_readsec("`cluster_var'", "`alluse'")
			mat `sec_vec_full' = r(sec_vec)
			
			** get unique values of `cluster_var'
			qui duplicates drop `cluster_var', force
			mata: s_readsec("`cluster_var'", "`alluse'")
			mat `sec_vec_unique' = r(sec_vec)
			
			** Compute SE
			mata: s_AKM0_cluster("`e_null'", "`ln'", "`Xdd'", "`Xddd'", "`Ydd'", `beta0', `critical_value', "`sec_vec_full'", "`sec_vec_unique'", `coef')
			
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

	** Save results for display
	ereturn scalar b = `coef'
	ereturn scalar se = `r(se)'
	ereturn scalar CI_upp = `r(CIu)'
	ereturn scalar CI_low = `r(CIl)'
	ereturn scalar p = `r(p)'
	ereturn scalar tstat = `r(tstat)'
	
end

*******************************************************************************
*************************** BEGIN MATA CODE ***********************************
*******************************************************************************

version 14.0
set matastrict on

mata:
void s_AKM1_consvec(string scalar regressor_var, 
					string scalar share_varlist, 
					string scalar dependant_var,
					string scalar touse)
{
	real matrix Mn
	real matrix ln
	real colvector tildeYn
	real matrix hat_theta
	real matrix hat_beta
	real scalar coef
	real matrix tildeZn
	real matrix tildeXn
	real matrix A
	real matrix Xdd
	real matrix Xddd
	real matrix Ydd
	real matrix e_ln
	
	// Generate Matrix of Regressors, shares, and outcome variable
	//Mn: Matrix of regressors
	//ln: Matrix of Shares
	//tildeYn: Vector of Dependent Variable
	
	st_view(Mn=., ., tokens(regressor_var), touse) 		
	st_view(ln=., ., tokens(share_varlist), touse) 		
	st_view(tildeYn=., ., tokens(dependant_var), touse) 
	
	hat_theta = invsym(Mn'*Mn) * (Mn' * tildeYn)
	e = tildeYn - Mn * hat_theta
	hat_beta = hat_theta[1, .]
	coef = hat_theta[1,1]
	
	tildeZn = Mn[|1,2\.,.|]
	tildeXn = Mn[|1,1\.,1|]
	
	A = tildeZn * invsym(tildeZn'*tildeZn)
	Ydd = tildeYn - A * (tildeZn'*tildeYn)
	Xdd = tildeXn - A * (tildeZn'*tildeXn)
	Xddd = invsym(ln' *ln) * (ln' * Xdd)
	e_ln = (e' * ln)'
	
	// convert from view to matrix that can be returned
	ln_q   = ln[|.,.\.,.|]
	Xdd_q  = Xdd[|.,.\.,.|]
	Xddd_q = Xddd[|.,.\.,.|]
	e_q 	 = e[|.,.\.,.|]
	e_ln_q 	 = e_ln[|.,.\.,.|]
	
	st_matrix("r(ln)", ln_q)
	st_matrix("r(Xdd)", Xdd_q)
	st_matrix("r(Xddd)", Xddd_q)
	st_matrix("r(e)", e_q)
	st_matrix("r(e_ln)", e_ln_q)
	st_numscalar("r(b)", coef)
}	

void s_AKM0_consvec(string scalar regressor_var, 
					string scalar share_varlist, 
					string scalar dependant_var,
					string scalar touse,
					scalar beta0)
{
	real matrix Mn
	real matrix ln
	real colvector tildeYn
	real matrix hat_theta
	real matrix hat_beta
	real scalar coef
	real matrix tildeZn
	real matrix tildeXn
	real matrix A
	real matrix Ydd
	real matrix Xdd
	real matrix Xddd
	real matrix e_ln
	real matrix e_null
	
	// Generate Matrix of Regressors, shares, and outcome variable
	//Mn: Matrix of regressors
	//ln: Matrix of Shares
	//tildeYn: Vector of Dependent Variable
	
	st_view(Mn=., ., tokens(regressor_var), touse) 		
	st_view(ln=., ., tokens(share_varlist), touse) 		
	st_view(tildeYn=., ., tokens(dependant_var), touse) 
	
	hat_theta = invsym(Mn'*Mn) * (Mn' * tildeYn)
	e = tildeYn - Mn * hat_theta
	hat_beta = hat_theta[1, .]
	coef = hat_theta[1,1]
	
	tildeZn = Mn[|1,2\.,.|]
	tildeXn = Mn[|1,1\.,1|]
	
	A = tildeZn * invsym(tildeZn'*tildeZn)
	Ydd = tildeYn - A * (tildeZn'*tildeYn)
	Xdd = tildeXn - A * (tildeZn'*tildeXn)
	Xddd = invsym(ln' *ln) * (ln' * Xdd)
	e_ln = (e' * ln)'
	e_null = Ydd - Xdd * beta0
	
	// convert from view to matrix that can be returned
	ln_q   	= ln[|.,.\.,.|]
	Xdd_q  	= Xdd[|.,.\.,.|]
	Xddd_q 	= Xddd[|.,.\.,.|]
	Ydd_q 	= Ydd[|.,.\.,.|]
	e_q 	= e[|.,.\.,.|]
	e_ln_q 	= e_ln[|.,.\.,.|]
	e_null_q = e_null[|.,.\.,.|]
	
	st_matrix("r(ln)", ln_q)
	st_matrix("r(Xdd)", Xdd_q)
	st_matrix("r(Xddd)", Xddd_q)
	st_matrix("r(Ydd)", Ydd_q)
	st_matrix("r(e)", e_q)
	st_matrix("r(e_ln)", e_ln_q)
	st_matrix("r(e_null)", e_null_q)
	st_numscalar("r(b)", coef)
}	
				  
void s_readsec(string scalar cluster_var, string scalar touse) 
{
	real matrix sec_vec
	
	st_view(sec_vec=., ., tokens(cluster_var), touse)
	
	sec_vec_q = sec_vec[|.,.\.,.|]
	
	st_matrix("r(sec_vec)", sec_vec_q)
}

void s_AKM1_nocluster(string scalar e_matrix,	 
					  string scalar ln_matrix,   
					  string scalar Xdd_matrix,  
					  string scalar Xddd_matrix, 
					  scalar coef, 				 
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
	variance = invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
	
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

void s_AKM0_nocluster(string scalar e_null_matrix,	 
					  string scalar ln_matrix,   
					  string scalar Xdd_matrix,  
					  string scalar Xddd_matrix, 
					  string scalar Ydd_matrix,
					  scalar coef, 		
					  scalar beta0,
					  scalar critical_value)
{	
	real matrix e_null
	real matrix ln
	real matrix Xdd
	real matrix Xddd
	real matrix Ydd
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
	
	// Co,pute variance matrix
	R = (e_null' * ln) :* (e_null' * ln)
	LambdaAKM = R * (Xddd :* Xddd)
	
	variance = invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
	
	se_AKMnull = sqrt(variance[1,1])
	tstat = (coef - beta0) / se_AKMnull
	
	p_value = 2 * (1 - normal(abs(tstat)))
	
	// Compute confidence interval
	critical2 = critical_value^2
    RY = Xdd' * Ydd
    RX = Xdd' * Xdd
	
	lnY =  Ydd' * ln 
    lnX =  Xdd' * ln
	
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
void s_AKM1_cluster(string scalar e_matrix,	 			
					string scalar ln_matrix,   			
					string scalar Xdd_matrix,  			
					string scalar Xddd_matrix, 			
					scalar beta0, 						
					scalar critical_value, 				
					string scalar sec_vec_full_matrix, 	
					string scalar sec_vec_unique_matrix, 
					scalar coef
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
	
	variance = invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
	
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
void s_AKM0_cluster(string scalar e_null_matrix,	 			
					string scalar ln_matrix,   			
					string scalar Xdd_matrix,  			
					string scalar Xddd_matrix, 	
					string scalar Ydd_matrix,
					scalar beta0, 						
					scalar critical_value, 				
					string scalar sec_vec_full_matrix, 	
					string scalar sec_vec_unique_matrix, 
					scalar coef
					)
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
        lnXCluster = (( Xdd'*lnCluster  )') :* XdddCluster
		
		LambdaAKM  = LambdaAKM  + sum(RXCluster * RXCluster')
		SXY = SXY + sum(lnYCluster * lnXCluster')
        SXX = SXX + sum(lnXCluster * lnXCluster')
        SYY = SYY + sum(lnYCluster * lnYCluster')
	}
	
	variance = invsym(Xdd' * Xdd) * LambdaAKM * invsym(Xdd' * Xdd)
	
	// AKM se, p-value, and t-stat
	se_AKM = sqrt(variance[1,1])
	tstat = (coef - beta0) / se_AKM
	
	p_value = 2 * (1 - normal(abs(tstat)))
	
	// Compute confidence interval
	critical2 = critical_value^2
    RY = Xdd' * Ydd
    RX = Xdd' * Xdd
	
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
