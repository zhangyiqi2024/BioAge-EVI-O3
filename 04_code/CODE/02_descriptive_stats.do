*****************************************************
* Descriptive Statistics Table (Baseline & Trajectory Groups, column headers with n)
* Output: ../output/tables/descriptive_stats_baseline.xlsx
*****************************************************
set seed 20240101

clear all
set more off

*—— 1) Load data and apply basic filters ——*
use "DATA/analysis_ready/analysis_dataset.dta", clear

keep if age > 45 & age < 90

*—— 2) Baseline sample only ——*
keep if year==2011

*—— 3) (Optional) Identify weight variable for weighted mean/SD of continuous vars ——*
local wvar ""
capture confirm variable weight
if _rc==0 local wvar "weight"
else {
    capture confirm variable weight
    if _rc==0 local wvar "weight"
    else {
        capture confirm variable pweight
        if _rc==0 local wvar "pweight"
    }
}

*—— 4) Variable list ——*
local convars  BAA_PhenoAge BA_PhenoAge EVI_1KM_mean4y O3_1KM_mean4y ///
               age TMP PER SNOW DEM lat bmi hh_consumption_total

local catvars  sex education retirement pension_insurance ///
               smoking_current alcohol_current hypertension dyslipidemia

*—— 5) Define groups (columns) ——*
gen byte g_all  = 1
gen byte g_evi0 = traj15_EVI_1KM==0  if !missing(traj15_EVI_1KM)
replace g_evi0  = 0 if missing(g_evi0)
gen byte g_evi1 = traj15_EVI_1KM==1  if !missing(traj15_EVI_1KM)
replace g_evi1  = 0 if missing(g_evi1)

gen byte g_o30  = traj15_O3_1KM==0   if !missing(traj15_O3_1KM)
replace g_o30   = 0 if missing(g_o30)
gen byte g_o31  = traj15_O3_1KM==1   if !missing(traj15_O3_1KM)
replace g_o31   = 0 if missing(g_o31)

*—— 5b) Count sample sizes per column (unweighted) ——*
quietly count if g_all==1
local N_all  = r(N)
quietly count if g_evi0==1
local N_evi0 = r(N)
quietly count if g_evi1==1
local N_evi1 = r(N)
quietly count if g_o30==1
local N_o30  = r(N)
quietly count if g_o31==1
local N_o31  = r(N)

*—— 6) Define cell-generation programs ——*
capture program drop _cell_cont_group
program define _cell_cont_group, rclass
    // Continuous variable: "Mean (SD)", can use pweight
    syntax varname, GROUP(varname) [WVAR(name)]
    tempvar use
    gen byte `use' = `group'==1 & !missing(`varlist')
    local wopt ""
    if "`wvar'"!="" local wopt "[pw=`wvar']"
    quietly summarize `varlist' `wopt' if `use'==1
    if (r(N)>0) {
        local m = string(round(r(mean), .01), "%9.2f")
        local s = string(round(r(sd),   .01), "%9.2f")
        return local cell = "`m' ("+"`s'"+")"
    }
    else    return local cell = ""
end

capture program drop _cell_cat_group_num
program define _cell_cat_group_num, rclass
    // Categorical variable (numeric): return "count (pct%)", denominator = group non-missing count (unweighted)
    syntax varname, GROUP(varname) LEVEL(real)
    tempvar use
    gen byte `use' = `group'==1
    quietly count if `use'==1 & !missing(`varlist')
    local denom = r(N)
    quietly count if `use'==1 & `varlist'==`level'
    local cnt = r(N)
    local pct = cond(`denom'>0, 100*`cnt'/`denom', .)
    if `denom'==0 return local cell = ""
    else return local cell = string(`cnt', "%9.0f") + " (" + string(round(`pct', .1), "%4.1f") + "%)"
end

capture program drop _cell_cat_group_str
program define _cell_cat_group_str, rclass
    // Categorical variable (string): return "count (pct%)", denominator = group non-missing count (unweighted)
    syntax varname, GROUP(varname) LEVEL(string)
    tempvar use
    gen byte `use' = `group'==1
    quietly count if `use'==1 & !missing(`varlist')
    local denom = r(N)
    quietly count if `use'==1 & `varlist'==`level'
    local cnt = r(N)
    local pct = cond(`denom'>0, 100*`cnt'/`denom', .)
    if `denom'==0 return local cell = ""
    else return local cell = string(`cnt', "%9.0f") + " (" + string(round(`pct', .1), "%4.1f") + "%)"
end

*—— 7) Begin collecting results ——*
tempfile T
tempname post
postfile `post' str80 rowname str60 all str60 evi0 str60 evi1 str60 o30 str60 o31 using `T', replace

* Continuous variables — one row per variable, 5 columns *
foreach v of local convars {
    * Use variable name directly as row label *
    local rn "`v'"

    foreach g in g_all g_evi0 g_evi1 g_o30 g_o31 {
        quietly _cell_cont_group `v', group(`g') wvar(`wvar')
        local c`g' "`r(cell)'"
    }
    post `post' ("`rn'") ("`cg_all'") ("`cg_evi0'") ("`cg_evi1'") ("`cg_o30'") ("`cg_o31'")
}

* Categorical variables — variable name row + level rows *
foreach v of local catvars {
    * Use variable name directly as row label *
    local vlab "`v'"

    * Variable name row (header only, data cells left blank)
    post `post' ("`vlab'") ("") ("") ("") ("") ("")

    * List all levels (clean to auto-quote strings with spaces)
    levelsof `v', local(levels) clean

    * Value label name (if present)
    local vl : value label `v'

    * For each level, append a row
    capture confirm numeric variable `v'
    local isnum = cond(_rc==0, 1, 0)

    foreach L of local levels {
        if `isnum' {
            local Ltxt "`L'"
            if "`vl'"!="" {
                local Ltxt : label `vl' `L'
            }
            local rn "  - `Ltxt'"
            foreach g in g_all g_evi0 g_evi1 g_o30 g_o31 {
                quietly _cell_cat_group_num `v', group(`g') level(`L')
                local c`g' "`r(cell)'"
            }
            post `post' ("`rn'") ("`cg_all'") ("`cg_evi0'") ("`cg_evi1'") ("`cg_o30'") ("`cg_o31'")
        }
        else {
            local Ltxt `"`L'"'
            local rn "  - `Ltxt'"
            foreach g in g_all g_evi0 g_evi1 g_o30 g_o31 {
                quietly _cell_cat_group_str `v', group(`g') level(`L')
                local c`g' "`r(cell)'"
            }
            post `post' ("`rn'") ("`cg_all'") ("`cg_evi0'") ("`cg_evi1'") ("`cg_o30'") ("`cg_o31'")
        }
    }
}

postclose `post'
use `T', clear

*—— 8) Insert column headers with n counts as first row ——*
insobs 1, before(1)
replace rowname = "Variable"                                     in 1
replace all     = "Baseline, all (n=`N_all')"                in 1
replace evi0    = "Low EVI trajectory (n=`N_evi0')"                   in 1
replace evi1    = "High EVI trajectory (n=`N_evi1')"                   in 1
replace o30     = "Low O3 trajectory (n=`N_o30')"                     in 1
replace o31     = "High O3 trajectory (n=`N_o31')"                     in 1

*—— 9) Export to Excel (without firstrow(var)) ——*
export excel using "output/tables/descriptive_stats_baseline.xlsx", sheet("Descriptive", replace)

display as result "Exported: ../output/tables/descriptive_stats_baseline.xlsx (headers include per-group sample size n)"
