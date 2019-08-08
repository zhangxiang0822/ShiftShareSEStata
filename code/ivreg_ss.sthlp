{smcl}
{* *! version 1.00  2aug2019}{...}
{viewerjumpto "Syntax" "ivreg_ss##syntax"}{...}
{viewerjumpto "Description" "ivreg_ss##description"}{...}
{viewerjumpto "Options" "ivreg_ss##options"}{...}
{viewerjumpto "Examples" "ivreg_ss##examples"}{...}
{viewerjumpto "Saved results" "ivreg_ss##saved_results"}{...}
{viewerjumpto "Author" "ivreg_ss##author"}{...}
{viewerjumpto "Author" "ivreg_ss##reference"}{...}
{title:Title}

{p2colset 5 18 24 2}{...}
{p2col :{hi: ivreg_ss} {hline 2}}Inference in IV regressions with shift-share structure. See also reg_ss for Inference in OLS regression with shift-share structure. 
{p2colreset}{...}

{marker syntax}{title:Syntax}

{p 8 15 2}
{cmd:ivreg_ss}
{depvar} 
{cmd:,} {opt endogenous_var(varlist)} {opt shiftshare_iv(varlist)} 
		{opt share_varlist(varlist)}  {opt alpha} {opt akmtype} 
		[control_varlist(varlist) weight_var(varlist) beta0 path_cluster cluster_var(varlist)] [{it:options}]

firststage(integer 0)]

{marker Options}{...}
{synoptset 29 tabbed}{...}
{synopthdr : Options}
{synoptline}
{syntab:Model}
{p2coldent:* {opt endogenous_var(varlist)}} Endogenous variable in the IV regression {p_end}
{p2coldent:* {opt shiftshare_iv(varlist)}} Shift-share variable {p_end}
{p2coldent:* {opt share_varlist(varlist)}} List of variable storing sector shares. One variable corresponds to one sector i's share in all regions. If there are s sectors, there should be s such variables. {p_end}
{p2coldent:* {opt alpha}} Determines confidence level of reported confidence intervals, which will have coverage 1-alpha. The default confidence level ia 5%. {p_end}
{p2coldent:* {opt akmtype}} Specifying which inference methods to use and can take the value 0 or 1. 1 for Adão-Kolesár-Morales (AKM) and 0 for AKM with null imposed. Note the reported standard
error for AKM0 method corresponds to the normalized standard error, given by the length of the confidence interval divided by 2*z(1-alpha/2). Also, if using AKM0 method,
you may also set beta0 value (default value for beta0 is 0) {p_end}

{synoptline}
{syntab:Other}
{synopt :{opt control_varlist(varlist)}} List of control variables to be included{p_end}
{synopt :{opt weight_var(varlist)}} An optional variable of weights to be used in the fitting process. If specified, for computing the first stage and the reduced
form, weighted least squares is used with weights "weights" (that is, minimizing sum(weights*residuals^2)); otherwise ordinary least squares is used. {p_end}
{synopt :{opt beta0 }} Null that is tested when using AKM0 method (only affects reported p-values) {p_end}
{synopt :{opt path_cluster}} The path to the .dta file which stores the cluster information (cluster-sector dataset). Provide this if you want to compute SE with sector clustered. {p_end}
{synopt :{opt cluster_var}} The variable name of the cluster variable in your cluster-sector dataset.{p_end}
{synopt :{opt firststage}}  Specifying whether to show first-stage results or not. It may take the value 1 if you want to see the first-stage results and inference, or 0 otherwise. The default value is 0. If you want to see first-stage results, you should also install reg_ss command{p_end}

{marker discription}{...}
{title:Description}
{pstd} {opt ivreg_ss} Provides standard error estimation and confidence intervals in instrumental variables regressions when the instrument has a shift-share structure.{p_end}

{marker examples}{...}
{title:Examples}
You may download the needed data at https://github.com/zhangxiang0822/BartikSEStata/tree/master/data

Example 1: AKM, no cluster
{phang2}{cmd:. 	use "data/ADH_derived.dta", clear }{p_end}
{phang2}{cmd:. 	local control_varlist t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource }{p_end}
{phang2}{cmd:.	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) control_varlist(`control_varlist') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) firststage(1)}{p_end}

Example 2: AKM0, with cluster
Be sure to put "sector_derived.dta" to your local path
{phang2}{cmd:. 	use "data/ADH_derived.dta", clear }{p_end}
{phang2}{cmd:. 	local control_varlist t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource }{p_end}
{phang2}{cmd:. 	local local_path "data/sector_derived.dta"}{p_end}
{phang2}{cmd:.	ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) control_varlist(`control_varlist') share_varlist(emp_share*) weight_var(weight) akmtype(0) path_cluster(`local_path') cluster_var(sec_3d)}

For more examples, please see "https://github.com/zhangxiang0822/BartikSEStata/blob/master/code/ADHapplication.do".

{marker author}{...}
{title:Bug Reporting}
{pstd} ivreg_ss is part of an ongoing project and thus may contain errors and malfunctions. Pleas submit bugs, comments, and suggestions to the Github repository at "https://github.com/zhangxiang0822/BartikSEStata". You're free to open an issue. {p_end}
			 
{marker author}{...}
{title:Author}
{pstd}Xiang Zhang{p_end}
{pstd}xiangzhang@princeton.edu{p_end}

{marker reference}{...}
{title:Reference}
{pstd} Adão, Rodrigo, Michal Kolesár, and Eduardo Morales. Shift-share designs: Theory and inference. No. w24944. National Bureau of Economic Research, 2018. {p_end}

{pstd} This code is adapted from MATLAB and R version code by Adão, Rodrigo, Michal Kolesár, and Eduardo Morales. {p_end}

{pstd} You may find MATLAB version code at "https://github.com/kolesarm/BartikSEMatlab". {p_end} 
{pstd} You may find R version code at "https://github.com/kolesarm/BartikSE". {p_end} 
