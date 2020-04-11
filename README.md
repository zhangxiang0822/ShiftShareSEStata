# ShiftShareSEStata
This Stata package implements confidence intervals in shift-share designs (also called [Bartik 1991](https://research.upjohn.org/up_press/77/) designs) using the `AKM` and `AKM0` procedures from [Adão, Kolesár, and Morales (2019)](https://doi.org/10.1093/qje/qjz025). See the
[ShiftShareSEMatlab](https://github.com/kolesarm/ShiftShareSEMatlab) package for
Matlab version of this code, and the
[ShiftShareSE](https://github.com/kolesarm/ShiftShareSE) package
for an R version.

## Example
IV regression using data from [Autor, Dorn, and Hanson (2013)](https://www.aeaweb.org/articles?id=10.1257/aer.103.6.2121), including the first-stage results

        use "data/ADH_derived.dta", clear
        local control_varlist t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource
        ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) control_varlist(`control_varlist') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) firststage(1)

Corresponding reduced-form regression

        reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) control_varlist(`control_varlist') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1)

The datasets needed to run these commands are stored in the `data` directory:
[`ADH_derived.dta`](data/ADH_derived.dta) and
[`sector_derived.dta`](data/sector_derived.dta). The do file
[`ADHapplication.do`](code/ADHapplication.do) contains additional examples that
use these datasets to illustrate how to use the `reg_ss` and `ivreg_ss` to
estimate linear regressions in which the regressor of interest has a shift-share
structure, and instrumental variables regressions in which the instrument has a
shift-share structure.

## Installation

You can install the released version of `reg_ss` and `ivreg_ss` from using the `ssc` command within Stata.

       ssc install reg_ss
       ssc install ivreg_ss

In the meantime, to install the package from Github, download the files
[`reg_ss.ado`](code/reg_ss.ado), [`reg_ss.sthlp`](code/reg_ss.sthlp),
[`ivreg_ss.ado`](code/ivreg_ss.ado), [`ivreg_ss.sthlp`](code/ivreg_ss.sthlp)
from the `code` directory, and put them into Stata's personal `ado` directory,
typically

- `c:\ado\personal` on Windows
- `~/Documents/Stata/ado` on a Mac
-  `~/ado` on Linux
- The simple way to find out your local Stata directory is to run the command `sysdir` in Stata.

For more information on how to use personal ado files, please refer to [Stata Official FAQ](https://www.stata.com/support/faqs/programming/personal-ado-directory/).

## Bug reporting and Questions
Please open issues and leave feedback, use click on the `Issues` tab.

## Some Special Notes
- Please make sure you have **no missing values** in your dataset. We don't
  handle missing values in our code, and having missing values in the dataset
  causes the program to shut down.
- Since Stata doesn't support nested `preserve` and `restore`, our `ivreg_ss`
  code cannot return your original dataset after computing AKM standard errors.
  `reg_ss` will return your original dataset if you don't compute clustered
  standard errors. We suggest reloading your dataset every time after you
  compute AKM standard errors.

## Common Error Message
- no observations: You may have missing values in your share variable list / shiftshare variable
