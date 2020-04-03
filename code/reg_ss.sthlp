{smcl}
{* *! version 1.00  2apr2020}{...}
{viewerjumpto "Syntax" "reg_ss##syntax"}{...}
{viewerjumpto "Options" "reg_ss##options"}{...}
{viewerjumpto "Collinear sectors and missing values" "reg_ss##description"}{...}
{viewerjumpto "Examples" "reg_ss##examples"}{...}
{viewerjumpto "Stored results" "reg_ss##results"}{...}
{viewerjumpto "Bug reporting" "reg_ss##bugs"}{...}
{viewerjumpto "Authors" "reg_ss##authors"}{...}
{viewerjumpto "References" "reg_ss##references"}{...}
{title:Title}

{p2colset 5 15 15 2}{...}
{p2col :{hi: reg_ss} {hline 2}}Computes confidence intervals, standard errors,
and p-values in a linear regression in which the regressor of interest has a
shift-share structure.
{p2colreset}{...}

{marker syntax}{title:Syntax}

{p 4 15 2}
{cmd:reg_ss} {depvar}{cmd:,} {opt shiftshare_var(varname)}
{opt share_varlist(varlist)}
[control_varlist(varlist) weight_var(varname) akmtype(#) beta0(#) path_cluster cluster_var(varname) alpha(#)]

{marker options}{...}
{synoptset 28}{...}
{synopthdr:Options}
{synoptline}
{syntab:Model}
{synopt :{opth shiftshare_var(varname)}}Shift-share regressor of interest.{p_end}
{synopt :{opth share_varlist(varlist)}}List of variables containing the sector
shares. Each variable indicates the regional shares for the corresponding
sector. Varlist must contain as many variables as there are sectors in the
analysis.{p_end}

{synoptline}
{syntab:Other options}
{synopt :{opt control_varlist(varlist)}}List of control variables to be included.{p_end}
{synopt :{opth weight_var(varname)}}Weights to be used in the fitting process.
If specified, weighted least squares is used;
otherwise, ordinary least squares is used.{p_end}
{synopt :{opt akmtype(#)}}Specifies which inference method to use. It must equal
0 or 1. Under {cmd:akmtype(1)}, {cmd:reg_ss} applies the AKM inference
procedure described in {help reg_ss##AKM:Adão, Kolesár, and Morales (2019)}.
Under {cmd:akmtype(0)}, {cmd:reg_ss} applies the AKM0 inference procedure
described in {help reg_ss##AKM:Adão, Kolesár, and Morales (2019)}. In this
case, the reported standard error corresponds to the normalized standard error,
given by the length of the corresponding confidence interval divided by
{cmd:2*invnormal(1-alpha/2)}. The default is {cmd:akmtype(1)}.{p_end}
{synopt :{opt beta0(#)}}Null that is tested when under {cmd:akmtype(0)} (only
affects reported p-values).{p_end}
{synopt :{opt path_cluster}}The path to the .dta file that stores the cluster
information corresponding to each sector. Provide this path if you want to
account for sectoral clusters.{p_end}
{synopt :{opth cluster_var(varname)}}The variable name of the cluster variable
in the .dta file storing the cluster information corresponding to each sector.
The order of the sectors in this .dta file should be the same as the order of
the sectors in the variable list containing the sector shares ({opt share_varlist(varlist)}).{p_end}
{synopt :{opt alpha(#)}}Determines confidence level of reported confidence
intervals, which will have coverage 1-{cmd:#}. The default value is
{cmd:alpha(0.05)}.{p_end}

{marker description}{...}
{title:Note on collinear share matrix}

{phang} If the share matrix implied by {cmd:share_varlist(varlist)} is
collinear, then it is not possible to to recover the sector-level shares, with
the controls partialled out (tilde{X}_s in the notation of
{help reg_ss##AKM:Adão, Kolesár, and Morales (2019)}), and {cmd:reg_ss} will
return an error message "Share matrix is collinear". To resolve this problem,
the researcher has three options:

{pmore}
(i) Drop the collinear sectors, and adjust the shift-share regressor (as
specified by {cmd:shiftshare_var(varname)}) accordinly. For example, suppose
that instead of the original S sectors, we only keep the first S0. Let w_is
denote the shares in location i and sector s, and let X_s denote the shock to
sector s, so that the original definition of the shift-share variable is X_i =
sum_{s=1}^{S}w_is*X_s. Then after we drop the collinear sectors, the new
shift-share variable will be given by X_i = sum_{s=1}^{S0}w_is*X_s. This
effectively puts shocks to the collinear sectors into the residual (which is
analogous to letting say the shock to non-manufacturing sectors be part of the
residual).

{pmore}
(ii) Aggregate the sectors. For example, instead of using 4-digit
manufacturing codes to define the sectors, use 3-digit codes, defining the shock
to the 3-digit sector as the sum of 4-digit sector shocks. Since the sectoral
shocks and shares are now different, we again need to adjust the definition of
shift-share regressor (as specified by
{cmd:shiftshare_var(varname)}).

{pmore}
(iii) If the only controls are those with shift-share structure, and we
have data on the sector-level variables Z_s that form these controls, we can
estimate tilde{X}_s directly by regressing the sectoral shocks on Z_s, and
taking the residual, and use this estimate in the standard error formula derived in
{help reg_ss##AKM:Adão, Kolesár, and Morales (2019)}. This third option is not
currently implemented in this package.
{p_end}

{phang} Options (i) and (ii) change the definition of the estimand, if there is
treatment effect heterogeneity. Since they involve changing the shift-share regressor
{cmd:shiftshare_var(varname)}), the adjustment to this variable needs to be done
before using the {cmd:reg_ss} command.
{p_end}

{title:Missing values}
{pstd} Before running regressions using {cmd:reg_ss},
please make sure there are no missing values in any variables that you pass to
the program.{p_end}

{marker examples}{title:Examples}
{phang}
The data for these examples may be downloaded from the
{browse "https://github.com/zhangxiang0822/ShiftShareSEStata/blob/master/data/ADH_derived.dta":GitHub repository} for this command. The data is a subset of data from {help reg_ss##ADH:Autor, Dorn, and Hanson (2013, ADH)}. The variables in this dataset are as follows:

{pmore}
{it:d_sh_empl}: Change in the share of working-age population.

{pmore}
{it:d_sh_empl_mfg}: Change in the share of working-age population employed in manufacturing.

{pmore}
{it:d_sh_empl_nmfg}: Change in the share of working-age population employed in non-manufacturing.

{pmore}
{it:d_tradeusch_pw}: Change in sectoral U.S. imports from China normalized by U.S. total employment in the corresponding sector, aggregated to regional level. This is the variable of interest in ADH.

{pmore}
{it:d_tradeotch_pw_lag}: Change in sectoral imports from China by by a set of high-income countries other than the U.S. normalized by the beginning-of-the period U.S. total employment in the corresponding sector, aggregated to the regional level. This is the instrumental variable in ADH.

{pmore}
{it:emp_share1-emp_share770}: The local employment share by region from the first to the 770th sector.

{pmore}
{it:weight}: Start-of-period share of U.S. population in each commuting zone (CZ).

{pmore}
{it:state}: State identification nunmber (FIPS code).

{pmore}
{it:czone}: CZ identification number.

{pmore}
{it:t2}: Indicator for the period 2000-2007, which corresponds to the second period in the sample.

{pmore}
{it:l_shind_manuf_cbp}: Manufacturing employment share in each CZ (in percentage terms).

{pmore}
{it:l_sh_popedu_c}: College-educated population share in each CZ (in percentage terms).

{pmore}
{it:l_sh_popfborn}: Foreign-born population share in each CZ (in percentage terms).

{pmore}
{it:l_sh_empl_f}: Female employment share in each CZ (in percentage terms).

{pmore}
{it:l_sh_routine33}: Share of employment in routine occupations in each CZ (in percentage terms).

{pmore}
{it:l_task_outsource}: Offshorability index of occupations in each CZ.
{p_end}

{pstd}Example 1: AKM standard errors without clustering{p_end}
{phang2}{cmd:.  use "data/ADH_derived.dta", clear }{p_end}
{phang2}{cmd:.  local controls t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource }{p_end}
{phang2}{cmd:.  reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) control_varlist(`controls') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1)}{p_end}

{pstd}Example 2: AKM0 standard errors without clustering.{p_end}
{phang2}{cmd:.  use "data/ADH_derived.dta", clear }{p_end}
{phang2}{cmd:.  local controls t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource }{p_end}
{phang2}{cmd:.  reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) control_varlist(`controls') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0)}

{pstd}Example 3: AKM standard errors, clustered at three digit SIC level. Make sure to put "sector_derived.dta" to your local path.{p_end}
{phang2}{cmd:.  use "data/ADH_derived.dta", clear }{p_end}
{phang2}{cmd:.  local controls t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource }{p_end}
{phang2}{cmd:.  local cpath "data/sector_derived.dta"}{p_end}
{phang2}{cmd:.  reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) control_varlist(`controls') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) path_cluster(`cpath') cluster_var(sec_3d) }

{pstd}Example 4: AKM0, clustered at three digit SIC level{p_end}
{phang2}{cmd:.  use "data/ADH_derived.dta", clear }{p_end}
{phang2}{cmd:.  local controls t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource }{p_end}
{phang2}{cmd:.  local cpath "data/sector_derived.dta"}{p_end}
{phang2}{cmd:.  reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) control_varlist(`controls') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(0) path_cluster(`cpath') cluster_var(sec_3d)}{p_end}

{pstd}More examples are available as a {browse "https://github.com/zhangxiang0822/BartikSEStata/blob/master/code/ADHapplication.do":do file} on the GitHub repository for this command.{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:reg_ss} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{synopt:{cmd:e(b)}}Estimate of the coefficient on the shift-share regressor {p_end}
{synopt:{cmd:e(se)}}Standard error of the estimate of the coefficient on the shift-share regressor {p_end}
{synopt:{cmd:e(CI_low)}}Lower bound of the confidence interval for the coefficient on the shift-share regressor{p_end}
{synopt:{cmd:e(CI_upp)}}Upper bound of the confidence interval for the coefficient on the shift-share regressor{p_end}
{synopt:{cmd:e(tstat)}}T-statistics for the null hypothesis that the coefficient on the shift-share regressor equals zero {p_end}
{synopt:{cmd:e(p)}}P-value for the null hypothesis that the coefficient on the shift-share regressor equals zero {p_end}

{marker bugs}{...}
{title:Bug Reporting}
{pstd}Please, submit bugs, comments, and suggestions by opening an issues at the GitHub repository at {browse "https://github.com/zhangxiang0822/BartikSEStata/issues"} for this command.{p_end}

{marker authors}{...}
{title:Authors}
{pstd}Rodrigo Adão{p_end}
{pstd}rodrigo.adao@chicagobooth.edu{p_end}

{pstd}Michal Kolesár{p_end}
{pstd}mkolesar@princeton.edu{p_end}

{pstd}Eduardo Morales{p_end}
{pstd}ecmorale@princeton.edu{p_end}

{pstd}Xiang Zhang{p_end}
{pstd}xiangzhang@princeton.edu{p_end}

{marker references}{...}
{title:References}

{marker AKM}{...}
{phang}
Adão, Rodrigo, Michal Kolesár, and Eduardo Morales. "Shift-share designs: Theory and inference." Quarterly Journal of Economics 134, no. 4 (2019): 1949-2010. {browse "https://doi.org/10.1093/qje/qjz025"}

{marker ADH}{...}
{phang}
Autor, David, David Dorn, and Gordon H. Hanson. "The China syndrome: Local labor market effects of import competition in the United States." American Economic Review 103, no. 6 (2013): 2121-2168. {browse "https://doi.org/10.1257/aer.103.6.2121"}

{phang}You may find a MATLAB version of this command at {browse "https://github.com/kolesarm/ShiftShareSEMatlab"}, and an R version at {browse "https://github.com/kolesarm/ShiftShareSE"}
{p_end}
