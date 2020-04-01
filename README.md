# ShiftShareSEStata
This R package implements confidence intervals in shift-share designs (also called [Bartik 1991](https://research.upjohn.org/up_press/77/) designs) using the `AKM` and `AKM0` procedures from [Adão, Kolesár, and Morales (2019)](https://doi.org/10.1093/qje/qjz025).

## Example
IV regression using data from [Autor, Dorn, and Hanson (2013)](https://www.aeaweb.org/articles?id=10.1257/aer.103.6.2121), including the first-stage results

        use "data/ADH_derived.dta", clear
        local control_varlist t2 l_shind_manuf_cbp reg_encen reg_escen reg_midatl reg_mount reg_pacif reg_satl reg_wncen reg_wscen l_sh_popedu_c l_sh_popfborn l_sh_empl_f l_sh_routine33 l_task_outsource
        ivreg_ss d_sh_empl_mfg, endogenous_var(d_tradeusch_pw) shiftshare_iv(d_tradeotch_pw_lag) control_varlist(`control_varlist') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1) firststage(1)

Corresponding reduced-form regression

        reg_ss d_sh_empl, shiftshare_var(d_tradeotch_pw_lag) control_varlist(`control_varlist') share_varlist(emp_share1-emp_share770) weight_var(weight) alpha(0.05) akmtype(1)

In `code/ADHapplication.do`, we provide additional examples on how to use the `reg_ss` and `ivreg_ss` to estimate linear regressions in which the regressor of interest/IV has a shift-share structure. To successfully run those codes, you may want to download the dataset needed from the `data` folder. The two datasets we use to produce results are `ADH_derived.dta` and `dector_derived.dta`.

## Installation
To install the package from Github, you may download the `reg_ss.ado`, `reg_ss.sthlp`, `ivreg_ss.ado`, and `ivreg_ss.sthlp`, and put them into the Stata personal ado directory.
- If you use Windows, it is probably `c:\ado\personal`, but it might be someplace else.
- If you use Mac, it is probably `~/Documents/Stata/ado`, but it might be someplace else.
- The simple way to find out your local Stata directory is to type `sysdir` in the Stata interface.

For more information on how to use personal ado files, please refer to [Stata Official FAQ](https://www.stata.com/support/faqs/programming/personal-ado-directory/).

In the near future, you will be able to directly install the packages from Stata.

## Bug reporting and Questions
You're welcomed to open issues and leave messages in the `Issues` section to report bugs or ask questions. 

## Some Special Notes
- Please make sure you have **no missing values** in your dataset. We don't handle missing-value problem in our code, and having missing values in the dataset would cause the probram shut down.
- Since Stata doesn't support nested `preserve` and `restore`, our `ivreg_ss` code cannot return your original dataset after computing AKM SE. `reg_ss` would return your original dataset if you don't compute clustered standard errors. We would suggest reload your dataset every time after you compute the AKM SE.

## Common Error Message
- no observations: You may have missing values in your share variable list/ shiftshare variable
