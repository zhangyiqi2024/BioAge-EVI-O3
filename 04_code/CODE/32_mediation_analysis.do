*==============================================================================
* Standardized multi-mediator analysis (Panel A/B/C + Nature-style tables)
* Pathway 1: log_EVI_1KM -> metabolic_activity -> BAA_PhenoAge
* Pathway 2: log_EVI_1KM -> O3_1KM -> BAA_PhenoAge
* Pathway 3: log_EVI_1KM -> social_participation -> BAA_PhenoAge
*==============================================================================

set seed 20240101

* Set up data
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90

* Variable transformations
gen ln_hh_consumption_total = log(hh_consumption_total)
gen log_EVI_1KM    = log(EVI_1KM_mean4y + 0.001)
gen log_EVI_500m   = log(EVI_500m_mean4y + 0.001)
gen log_O3_1KM     = log(O3_1KM_mean4y + 0.001)
gen log_O3_500m    = log(O3_500m_mean4y + 0.001)

* Principal component analysis for climatic variables
di "=== Principal component analysis for climatic variables ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* Define control variable global macro
global controls "c.age i.sex i.education c.ln_hh_consumption_total i.retirement i.pension_insurance i.smoking_current i.alcohol_current climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 c.bmi i.hypertension i.dyslipidemia i.year"

* ==============================================================================
* Data standardization
* ==============================================================================

di " "
di "=== Standardization of key variables ==="

* Standardize key variables (with protection for sd=0)
foreach var in BAA_PhenoAge metabolic_activity O3_1KM social_participation log_EVI_1KM {
    qui sum `var'
    if r(sd) > 0 {
        gen `var'_std = (`var' - r(mean)) / r(sd)
    }
    else {
        gen `var'_std = .
        di as error "Warning: variable `var' has zero standard deviation; cannot be standardized."
    }
    di "Variable `var' standardized: mean = " r(mean) ", SD = " r(sd)
}

* ==============================================================================
* Standardized mediation analysis 1: log_EVI_1KM -> metabolic_activity -> BAA_PhenoAge
* ==============================================================================

di " "
di "**********************************************************************"
di "Standardized mediation analysis 1: log_EVI_1KM -> total metabolic activity -> BAA_PhenoAge"
di "**********************************************************************"

* Clear any existing scalars
capture scalar drop a1_std b1_std c_total1_std c_direct1_std indirect1_std ///
    se_a1_std se_b1_std se_c_total1_std se_c_direct1_std se_indirect1_std ///
    p_a1_std p_b1_std p_c_total1_std p_c_direct1_std p_indirect1_std ///
    z_indirect1_std prop_mediated1_std z_c_total1_std z_c_direct1_std

* Step 1: total effect (X -> Y), using sample with metabolic_activity data
di "=== Step 1: total effect (X -> Y) for metabolic activity pathway ==="

mixed BAA_PhenoAge_std c.log_EVI_1KM_std $controls ///
    if !missing(metabolic_activity_std) ///
    || id_num:, vce(robust)

scalar c_total1_std    = _b[log_EVI_1KM_std]
scalar se_c_total1_std = _se[log_EVI_1KM_std]
scalar z_c_total1_std  = c_total1_std / se_c_total1_std
scalar p_c_total1_std  = 2*(1 - normal(abs(z_c_total1_std)))

est store step1_metabolic_std

di "Standardized total effect c (metabolic pathway) = " ///
   %6.3f c_total1_std ", SE = " %6.3f se_c_total1_std ", p-value = " %6.3f p_c_total1_std

* Step 2: X -> M (exposure to mediator)
di "=== Step 2: path a (X -> mediator) ==="

mixed metabolic_activity_std c.log_EVI_1KM_std $controls ///
     || id_num:, vce(robust)

scalar a1_std    = _b[log_EVI_1KM_std]
scalar se_a1_std = _se[log_EVI_1KM_std]
scalar z_a1_std  = a1_std / se_a1_std
scalar p_a1_std  = 2*(1 - normal(abs(z_a1_std)))

est store step2_metabolic_std

di "Standardized path a = " %6.3f a1_std ///
   ", SE = " %6.3f se_a1_std ", p-value = " %6.3f p_a1_std

* Step 3: direct effect (X and mediator jointly in the model)
di "=== Step 3: direct effect (X and mediator jointly in the model) ==="

mixed BAA_PhenoAge_std c.log_EVI_1KM_std c.metabolic_activity_std $controls ///
     || id_num:, vce(robust)

scalar b1_std    = _b[metabolic_activity_std]
scalar se_b1_std = _se[metabolic_activity_std]
scalar z_b1_std  = b1_std / se_b1_std
scalar p_b1_std  = 2*(1 - normal(abs(z_b1_std)))

scalar c_direct1_std    = _b[log_EVI_1KM_std]
scalar se_c_direct1_std = _se[log_EVI_1KM_std]
scalar z_c_direct1_std  = c_direct1_std / se_c_direct1_std
scalar p_c_direct1_std  = 2*(1 - normal(abs(z_c_direct1_std)))

est store step3_metabolic_std

di "Standardized path b = " %6.3f b1_std ///
   ", SE = " %6.3f se_b1_std ", p-value = " %6.3f p_b1_std
di "Standardized direct effect c' = " %6.3f c_direct1_std ///
   ", SE = " %6.3f se_c_direct1_std ", p-value = " %6.3f p_c_direct1_std

* Calculate standardized indirect effect
di "=== Indirect effect (standardized) ==="
scalar indirect1_std    = a1_std * b1_std
scalar se_indirect1_std = sqrt(a1_std^2 * se_b1_std^2 + b1_std^2 * se_a1_std^2)
scalar z_indirect1_std  = indirect1_std / se_indirect1_std
scalar p_indirect1_std  = 2 * (1 - normal(abs(z_indirect1_std)))

if c_total1_std != 0 {
    scalar prop_mediated1_std = indirect1_std / c_total1_std
}
else {
    scalar prop_mediated1_std = .
    di "Warning: total effect is zero; proportion mediated cannot be computed."
}

* ==============================================================================
* Standardized mediation analysis 2: log_EVI_1KM -> O3_1KM -> BAA_PhenoAge
* ==============================================================================

di " "
di "**********************************************************************"
di "Standardized mediation analysis 2: log_EVI_1KM -> O3_1KM -> BAA_PhenoAge"
di "**********************************************************************"

* Clear any existing scalars
capture scalar drop a2_std b2_std c_total2_std c_direct2_std indirect2_std ///
    se_a2_std se_b2_std se_c_total2_std se_c_direct2_std se_indirect2_std ///
    p_a2_std p_b2_std p_c_total2_std p_c_direct2_std p_indirect2_std ///
    z_indirect2_std prop_mediated2_std z_c_total2_std z_c_direct2_std

* Step 1: total effect (X -> Y), using sample with O3 data
di "=== Step 1: total effect (X -> Y) for O3 pathway ==="

mixed BAA_PhenoAge_std c.log_EVI_1KM_std $controls ///
    if !missing(O3_1KM_std) ///
    || id_num:, vce(robust)

scalar c_total2_std    = _b[log_EVI_1KM_std]
scalar se_c_total2_std = _se[log_EVI_1KM_std]
scalar z_c_total2_std  = c_total2_std / se_c_total2_std
scalar p_c_total2_std  = 2*(1 - normal(abs(z_c_total2_std)))

est store step1_pollution_std

di "Standardized total effect c (O3 pathway) = " ///
   %6.3f c_total2_std ", SE = " %6.3f se_c_total2_std ", p-value = " %6.3f p_c_total2_std

* Step 2: X -> M (exposure to mediator)
di "=== Step 2: path a (X -> mediator) ==="

mixed O3_1KM_std c.log_EVI_1KM_std $controls ///
     || id_num:, vce(robust)

scalar a2_std    = _b[log_EVI_1KM_std]
scalar se_a2_std = _se[log_EVI_1KM_std]
scalar z_a2_std  = a2_std / se_a2_std
scalar p_a2_std  = 2*(1 - normal(abs(z_a2_std)))

est store step2_pollution_std

di "Standardized path a = " %6.3f a2_std ///
   ", SE = " %6.3f se_a2_std ", p-value = " %6.3f p_a2_std

* Step 3: direct effect (X and mediator jointly in the model)
di "=== Step 3: direct effect (X and mediator jointly in the model) ==="

mixed BAA_PhenoAge_std c.log_EVI_1KM_std c.O3_1KM_std $controls ///
     || id_num:, vce(robust)

scalar b2_std    = _b[O3_1KM_std]
scalar se_b2_std = _se[O3_1KM_std]
scalar z_b2_std  = b2_std / se_b2_std
scalar p_b2_std  = 2*(1 - normal(abs(z_b2_std)))

scalar c_direct2_std    = _b[log_EVI_1KM_std]
scalar se_c_direct2_std = _se[log_EVI_1KM_std]
scalar z_c_direct2_std  = c_direct2_std / se_c_direct2_std
scalar p_c_direct2_std  = 2*(1 - normal(abs(z_c_direct2_std)))

est store step3_pollution_std

di "Standardized path b = " %6.3f b2_std ///
   ", SE = " %6.3f se_b2_std ", p-value = " %6.3f p_b2_std
di "Standardized direct effect c' = " %6.3f c_direct2_std ///
   ", SE = " %6.3f se_c_direct2_std ", p-value = " %6.3f p_c_direct2_std

* Calculate standardized indirect effect
di "=== Indirect effect (standardized) ==="
scalar indirect2_std    = a2_std * b2_std
scalar se_indirect2_std = sqrt(a2_std^2 * se_b2_std^2 + b2_std^2 * se_a2_std^2)
scalar z_indirect2_std  = indirect2_std / se_indirect2_std
scalar p_indirect2_std  = 2 * (1 - normal(abs(z_indirect2_std)))

if c_total2_std != 0 {
    scalar prop_mediated2_std = indirect2_std / c_total2_std
}
else {
    scalar prop_mediated2_std = .
}

* ==============================================================================
* Standardized mediation analysis 3: log_EVI_1KM -> social_participation -> BAA_PhenoAge
* ==============================================================================

di " "
di "**********************************************************************"
di "Standardized mediation analysis 3: log_EVI_1KM -> social participation -> BAA_PhenoAge"
di "**********************************************************************"

* Clear any existing scalars
capture scalar drop a3_std b3_std c_total3_std c_direct3_std indirect3_std ///
    se_a3_std se_b3_std se_c_total3_std se_c_direct3_std se_indirect3_std ///
    p_a3_std p_b3_std p_c_total3_std p_c_direct3_std p_indirect3_std ///
    z_indirect3_std prop_mediated3_std z_c_total3_std z_c_direct3_std

* Step 1: total effect (X -> Y), using sample with social participation data
di "=== Step 1: total effect (X -> Y) for social participation pathway ==="

mixed BAA_PhenoAge_std c.log_EVI_1KM_std $controls ///
    if !missing(social_participation_std) ///
    || id_num:, vce(robust)

scalar c_total3_std    = _b[log_EVI_1KM_std]
scalar se_c_total3_std = _se[log_EVI_1KM_std]
scalar z_c_total3_std  = c_total3_std / se_c_total3_std
scalar p_c_total3_std  = 2*(1 - normal(abs(z_c_total3_std)))

est store step1_social_std

di "Standardized total effect c (social participation pathway) = " ///
   %6.3f c_total3_std ", SE = " %6.3f se_c_total3_std ", p-value = " %6.3f p_c_total3_std

* Step 2: X -> M (exposure to mediator)
di "=== Step 2: path a (X -> mediator) ==="

mixed social_participation_std c.log_EVI_1KM_std $controls ///
     || id_num:, vce(robust)

scalar a3_std    = _b[log_EVI_1KM_std]
scalar se_a3_std = _se[log_EVI_1KM_std]
scalar z_a3_std  = a3_std / se_a3_std
scalar p_a3_std  = 2*(1 - normal(abs(z_a3_std)))

est store step2_social_std

di "Standardized path a = " %6.3f a3_std ///
   ", SE = " %6.3f se_a3_std ", p-value = " %6.3f p_a3_std

* Step 3: direct effect (X and mediator jointly in the model)
di "=== Step 3: direct effect (X and mediator jointly in the model) ==="

mixed BAA_PhenoAge_std c.log_EVI_1KM_std c.social_participation_std $controls ///
     || id_num:, vce(robust)

scalar b3_std    = _b[social_participation_std]
scalar se_b3_std = _se[social_participation_std]
scalar z_b3_std  = b3_std / se_b3_std
scalar p_b3_std  = 2*(1 - normal(abs(z_b3_std)))

scalar c_direct3_std    = _b[log_EVI_1KM_std]
scalar se_c_direct3_std = _se[log_EVI_1KM_std]
scalar z_c_direct3_std  = c_direct3_std / se_c_direct3_std
scalar p_c_direct3_std  = 2*(1 - normal(abs(z_c_direct3_std)))

est store step3_social_std

di "Standardized path b = " %6.3f b3_std ///
   ", SE = " %6.3f se_b3_std ", p-value = " %6.3f p_b3_std
di "Standardized direct effect c' = " %6.3f c_direct3_std ///
   ", SE = " %6.3f se_c_direct3_std ", p-value = " %6.3f p_c_direct3_std

* Calculate standardized indirect effect
di "=== Indirect effect (standardized) ==="
scalar indirect3_std    = a3_std * b3_std
scalar se_indirect3_std = sqrt(a3_std^2 * se_b3_std^2 + b3_std^2 * se_a3_std^2)
scalar z_indirect3_std  = indirect3_std / se_indirect3_std
scalar p_indirect3_std  = 2 * (1 - normal(abs(z_indirect3_std)))

if c_total3_std != 0 {
    scalar prop_mediated3_std = indirect3_std / c_total3_std
}
else {
    scalar prop_mediated3_std = .
}

* ==============================================================================
* Brief summary (for terminal display)
* ==============================================================================

di " "
di "**********************************************************************"
di "Standardized multi-mediator analysis: summary of path coefficients"
di "**********************************************************************"

di " "
di "Table 1. Standardized multi-mediator analysis: path coefficients"
di "========================================================================================================================"
di "Mediating pathway             X->M coefficient (a)   M->Y coefficient (b)   Total effect (c)   Indirect effect (axb)   Proportion (%)   Sobel z   p-value"
di "------------------------------------------------------------------------------------------------------------------------"

foreach i in 1 2 3 {

    capture scalar list a`i'_std
    if _rc == 0 {
        local a_val   = a`i'_std
        local b_val   = b`i'_std
        local c_val   = c_total`i'_std
        local ind_val = indirect`i'_std
        local prop_val = prop_mediated`i'_std * 100
        local z_val   = z_indirect`i'_std
        local p_val   = p_indirect`i'_std

        local stars_a   = cond(p_a`i'_std<0.001,"***",cond(p_a`i'_std<0.01,"**",cond(p_a`i'_std<0.05,"*","")))
        local stars_b   = cond(p_b`i'_std<0.001,"***",cond(p_b`i'_std<0.01,"**",cond(p_b`i'_std<0.05,"*","")))
        local stars_c   = cond(p_c_total`i'_std<0.001,"***",cond(p_c_total`i'_std<0.01,"**",cond(p_c_total`i'_std<0.05,"*","")))
        local stars_ind = cond(p_indirect`i'_std<0.001,"***",cond(p_indirect`i'_std<0.01,"**",cond(p_indirect`i'_std<0.05,"*","")))

        if `i' == 1 {
            local path_name "Metabolic activity pathway"
        }
        else if `i' == 2 {
            local path_name "Ambient ozone exposure pathway"
        }
        else if `i' == 3 {
            local path_name "Social participation pathway"
        }

        di "`path_name'" _col(30) %10.4f `a_val' "`stars_a'" _col(53) ///
           %10.4f `b_val' "`stars_b'" _col(76) ///
           %10.4f `c_val' "`stars_c'" _col(99) ///
           %10.4f `ind_val' _col(122) ///
           %6.2f `prop_val' _col(137) ///
           %10.4f `z_val' _col(152) ///
           %10.4f `p_val' "`stars_ind'"
    }
}

di "========================================================================================================================"
di "Note: *** p<0.001, ** p<0.01, * p<0.05; all coefficients are standardized."

* Output direct effect information for each pathway
di " "
di "Table 2. Standardized direct effects (c')"
di "=================================================================="
di "Mediating pathway             Direct effect (c')   SE        p-value"
di "------------------------------------------------------------------"

foreach i in 1 2 3 {
    capture scalar list c_direct`i'_std
    if _rc == 0 {
        local dir_val = c_direct`i'_std
        local se_val  = se_c_direct`i'_std
        local p_val   = p_c_direct`i'_std

        local stars_dir = cond(p_c_direct`i'_std<0.001,"***",cond(p_c_direct`i'_std<0.01,"**",cond(p_c_direct`i'_std<0.05,"*","")))

        if `i' == 1 {
            local path_name "Metabolic activity pathway"
        }
        else if `i' == 2 {
            local path_name "Ambient ozone exposure pathway"
        }
        else if `i' == 3 {
            local path_name "Social participation pathway"
        }

        di "`path_name'" _col(30) %10.4f `dir_val' _col(47) ///
           %8.4f `se_val' _col(60) ///
           %8.4f `p_val' "`stars_dir'"
    }
}

di "=================================================================="

* ==============================================================================
* Save log file (wide format, for internal archival)
* ==============================================================================

di "=== Save results to text log (wide format, internal use only) ==="
capture log close multi_mediation_results_std
log using "output/tables/mediation_results_summary.txt", text replace name(multi_mediation_results_std)

di "Standardized multi-mediator analysis summary"
di "Analysis time: $S_TIME $S_DATE"
di "Exposure: log_EVI_1KM; Outcome: PhenoAge-based biological age acceleration (BAA_PhenoAge)"
di "Sample size (depends on model and mediator): see mixed-model output."
di " "

di "Table 1. Standardized multi-mediator analysis: path coefficients"
di "========================================================================================================================"
di "Mediating pathway             X->M coefficient (a)   M->Y coefficient (b)   Total effect (c)   Indirect effect (axb)   Proportion (%)   Sobel z   p-value"
di "------------------------------------------------------------------------------------------------------------------------"

foreach i in 1 2 3 {
    capture scalar list a`i'_std
    if _rc == 0 {
        local a_val   = a`i'_std
        local b_val   = b`i'_std
        local c_val   = c_total`i'_std
        local ind_val = indirect`i'_std
        local prop_val = prop_mediated`i'_std * 100
        local z_val   = z_indirect`i'_std
        local p_val   = p_indirect`i'_std

        local stars_a   = cond(p_a`i'_std<0.001,"***",cond(p_a`i'_std<0.01,"**",cond(p_a`i'_std<0.05,"*","")))
        local stars_b   = cond(p_b`i'_std<0.001,"***",cond(p_b`i'_std<0.01,"**",cond(p_b`i'_std<0.05,"*","")))
        local stars_c   = cond(p_c_total`i'_std<0.001,"***",cond(p_c_total`i'_std<0.01,"**",cond(p_c_total`i'_std<0.05,"*","")))
        local stars_ind = cond(p_indirect`i'_std<0.001,"***",cond(p_indirect`i'_std<0.01,"**",cond(p_indirect`i'_std<0.05,"*","")))

        if `i' == 1 {
            local path_name "Metabolic activity pathway"
        }
        else if `i' == 2 {
            local path_name "Ambient ozone exposure pathway"
        }
        else if `i' == 3 {
            local path_name "Social participation pathway"
        }

        di "`path_name'" _col(30) %10.4f `a_val' "`stars_a'" _col(53) ///
           %10.4f `b_val' "`stars_b'" _col(76) ///
           %10.4f `c_val' "`stars_c'" _col(99) ///
           %10.4f `ind_val' _col(122) ///
           %6.2f `prop_val' _col(137) ///
           %10.4f `z_val' _col(152) ///
           %10.4f `p_val' "`stars_ind'"
    }
}

di "========================================================================================================================"
di "Note: *** p<0.001, ** p<0.01, * p<0.05; all coefficients are standardized."

log close multi_mediation_results_std

* Save estimation results
estimates save "output/tables/mediation_results", replace

di " "
di "=== Estimation results saved ==="
di "Results file saved: ../output/tables/mediation_results_summary.txt"
di "Estimation results saved: ../output/tables/mediation_results.ster"

* ==============================================================================
* Panel A/B/C panel mediation results table (Total effect (c) reported only in Panel A)
*   Output file: mediation_analysis_results.csv
* ==============================================================================

di " "
di "=== Generating compact Panel A/B/C mediation results table (long format) ==="

preserve
clear
set obs 16   // Panel A: 6 rows; Panel B: 5 rows; Panel C: 5 rows (16 rows total)

gen str35 panel   = ""
gen str30 effect  = ""
gen str30 beta_se = ""
gen str30 ci95    = ""
gen str10 p_value = ""
gen str10 pm_pct  = ""

local row = 0

*-------------------------------
* Panel A: Metabolic activity
*-------------------------------
local panel_name "Panel A: Metabolic activity"

* a: X->M
local row = `row' + 1
replace panel  = "`panel_name'" in `row'
replace effect = "a: X->M"        in `row'
local ba  = a1_std
local sea = se_a1_std
local loa = `ba' - 1.96*`sea'
local hia = `ba' + 1.96*`sea'
local pa  = p_a1_std
replace beta_se = string(`ba',"%6.3f") + " (" + string(`sea',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`loa',"%6.3f") + ", " + string(`hia',"%6.3f") + "]" in `row'
if `pa' < 0.001      replace p_value = "<0.001" in `row'
else                 replace p_value = string(`pa',"%6.3f") in `row'

* b: M->Y
local row = `row' + 1
replace panel  = "`panel_name'" in `row'
replace effect = "b: M->Y"        in `row'
local bb  = b1_std
local seb = se_b1_std
local lob = `bb' - 1.96*`seb'
local hib = `bb' + 1.96*`seb'
local pb  = p_b1_std
replace beta_se = string(`bb',"%6.3f") + " (" + string(`seb',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`lob',"%6.3f") + ", " + string(`hib',"%6.3f") + "]" in `row'
if `pb' < 0.001      replace p_value = "<0.001" in `row'
else                 replace p_value = string(`pb',"%6.3f") in `row'

* Total effect (c) -- reported only in Panel A
local row = `row' + 1
replace panel  = "`panel_name'"     in `row'
replace effect = "Total effect (c)" in `row'
local ct   = c_total1_std
local sect = se_c_total1_std
local loct = `ct' - 1.96*`sect'
local hict = `ct' + 1.96*`sect'
local pct  = p_c_total1_std
replace beta_se = string(`ct',"%6.3f") + " (" + string(`sect',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`loct',"%6.3f") + ", " + string(`hict',"%6.3f") + "]" in `row'
if `pct' < 0.001      replace p_value = "<0.001" in `row'
else                  replace p_value = string(`pct',"%6.3f") in `row'

* Direct effect (c')
local row = `row' + 1
replace panel  = "`panel_name'"       in `row'
replace effect = "Direct effect (c')" in `row'
local cd   = c_direct1_std
local secd = se_c_direct1_std
local locd = `cd' - 1.96*`secd'
local hicd = `cd' + 1.96*`secd'
local pcd  = p_c_direct1_std
replace beta_se = string(`cd',"%6.3f") + " (" + string(`secd',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`locd',"%6.3f") + ", " + string(`hicd',"%6.3f") + "]" in `row'
if `pcd' < 0.001      replace p_value = "<0.001" in `row'
else                  replace p_value = string(`pcd',"%6.3f") in `row'

* Indirect effect (axb)
local row = `row' + 1
replace panel  = "`panel_name'"           in `row'
replace effect = "Indirect effect (axb)"  in `row'
local ind   = indirect1_std
local seind = se_indirect1_std
local loind = `ind' - 1.96*`seind'
local hiind = `ind' + 1.96*`seind'
local pind  = p_indirect1_std
replace beta_se = string(`ind',"%6.3f") + " (" + string(`seind',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`loind',"%6.3f") + ", " + string(`hiind',"%6.3f") + "]" in `row'
if `pind' < 0.001      replace p_value = "<0.001" in `row'
else                   replace p_value = string(`pind',"%6.3f") in `row'

* Proportion mediated (%)
local row = `row' + 1
replace panel  = "`panel_name'"              in `row'
replace effect = "Proportion mediated (%)"   in `row'
local pm1 = prop_mediated1_std*100
if !missing(`pm1') replace pm_pct = string(`pm1',"%6.3f") in `row'

*-------------------------------
* Panel B: Ambient ozone exposure (no repeat of Total effect (c))
*-------------------------------
local panel_name "Panel B: Ambient ozone exposure"

* a: X->M
local row = `row' + 1
replace panel  = "`panel_name'" in `row'
replace effect = "a: X->M"        in `row'
local ba2  = a2_std
local sea2 = se_a2_std
local loa2 = `ba2' - 1.96*`sea2'
local hia2 = `ba2' + 1.96*`sea2'
local pa2  = p_a2_std
replace beta_se = string(`ba2',"%6.3f") + " (" + string(`sea2',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`loa2',"%6.3f") + ", " + string(`hia2',"%6.3f") + "]" in `row'
if `pa2' < 0.001      replace p_value = "<0.001" in `row'
else                  replace p_value = string(`pa2',"%6.3f") in `row'

* b: M->Y
local row = `row' + 1
replace panel  = "`panel_name'" in `row'
replace effect = "b: M->Y"        in `row'
local bb2  = b2_std
local seb2 = se_b2_std
local lob2 = `bb2' - 1.96*`seb2'
local hib2 = `bb2' + 1.96*`seb2'
local pb2  = p_b2_std
replace beta_se = string(`bb2',"%6.3f") + " (" + string(`seb2',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`lob2',"%6.3f") + ", " + string(`hib2',"%6.3f") + "]" in `row'
if `pb2' < 0.001      replace p_value = "<0.001" in `row'
else                  replace p_value = string(`pb2',"%6.3f") in `row'

* Direct effect (c')
local row = `row' + 1
replace panel  = "`panel_name'"       in `row'
replace effect = "Direct effect (c')" in `row'
local cd2   = c_direct2_std
local secd2 = se_c_direct2_std
local locd2 = `cd2' - 1.96*`secd2'
local hicd2 = `cd2' + 1.96*`secd2'
local pcd2  = p_c_direct2_std
replace beta_se = string(`cd2',"%6.3f") + " (" + string(`secd2',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`locd2',"%6.3f") + ", " + string(`hicd2',"%6.3f") + "]" in `row'
if `pcd2' < 0.001      replace p_value = "<0.001" in `row'
else                   replace p_value = string(`pcd2',"%6.3f") in `row'

* Indirect effect (axb)
local row = `row' + 1
replace panel  = "`panel_name'"           in `row'
replace effect = "Indirect effect (axb)"  in `row'
local ind2   = indirect2_std
local seind2 = se_indirect2_std
local loind2 = `ind2' - 1.96*`seind2'
local hiind2 = `ind2' + 1.96*`seind2'
local pind2  = p_indirect2_std
replace beta_se = string(`ind2',"%6.3f") + " (" + string(`seind2',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`loind2',"%6.3f") + ", " + string(`hiind2',"%6.3f") + "]" in `row'
if `pind2' < 0.001      replace p_value = "<0.001" in `row'
else                    replace p_value = string(`pind2',"%6.3f") in `row'

* Proportion mediated (%)
local row = `row' + 1
replace panel  = "`panel_name'"              in `row'
replace effect = "Proportion mediated (%)"   in `row'
local pm2 = prop_mediated2_std*100
if !missing(`pm2') replace pm_pct = string(`pm2',"%6.3f") in `row'

*-------------------------------
* Panel C: Social participation (likewise no repeat of c)
*-------------------------------
local panel_name "Panel C: Social participation"

* a: X->M
local row = `row' + 1
replace panel  = "`panel_name'" in `row'
replace effect = "a: X->M"        in `row'
local ba3  = a3_std
local sea3 = se_a3_std
local loa3 = `ba3' - 1.96*`sea3'
local hia3 = `ba3' + 1.96*`sea3'
local pa3  = p_a3_std
replace beta_se = string(`ba3',"%6.3f") + " (" + string(`sea3',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`loa3',"%6.3f") + ", " + string(`hia3',"%6.3f") + "]" in `row'
if `pa3' < 0.001      replace p_value = "<0.001" in `row'
else                  replace p_value = string(`pa3',"%6.3f") in `row'

* b: M->Y
local row = `row' + 1
replace panel  = "`panel_name'" in `row'
replace effect = "b: M->Y"        in `row'
local bb3  = b3_std
local seb3 = se_b3_std
local lob3 = `bb3' - 1.96*`seb3'
local hib3 = `bb3' + 1.96*`seb3'
local pb3  = p_b3_std
replace beta_se = string(`bb3',"%6.3f") + " (" + string(`seb3',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`lob3',"%6.3f") + ", " + string(`hib3',"%6.3f") + "]" in `row'
if `pb3' < 0.001      replace p_value = "<0.001" in `row'
else                  replace p_value = string(`pb3',"%6.3f") in `row'

* Direct effect (c')
local row = `row' + 1
replace panel  = "`panel_name'"       in `row'
replace effect = "Direct effect (c')" in `row'
local cd3   = c_direct3_std
local secd3 = se_c_direct3_std
local locd3 = `cd3' - 1.96*`secd3'
local hicd3 = `cd3' + 1.96*`secd3'
local pcd3  = p_c_direct3_std
replace beta_se = string(`cd3',"%6.3f") + " (" + string(`secd3',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`locd3',"%6.3f") + ", " + string(`hicd3',"%6.3f") + "]" in `row'
if `pcd3' < 0.001      replace p_value = "<0.001" in `row'
else                   replace p_value = string(`pcd3',"%6.3f") in `row'

* Indirect effect (axb)
local row = `row' + 1
replace panel  = "`panel_name'"           in `row'
replace effect = "Indirect effect (axb)"  in `row'
local ind3   = indirect3_std
local seind3 = se_indirect3_std
local loind3 = `ind3' - 1.96*`seind3'
local hiind3 = `ind3' + 1.96*`seind3'
local pind3  = p_indirect3_std
replace beta_se = string(`ind3',"%6.3f") + " (" + string(`seind3',"%6.3f") + ")" in `row'
replace ci95    = "[" + string(`loind3',"%6.3f") + ", " + string(`hiind3',"%6.3f") + "]" in `row'
if `pind3' < 0.001      replace p_value = "<0.001" in `row'
else                    replace p_value = string(`pind3',"%6.3f") in `row'

* Proportion mediated (%)
local row = `row' + 1
replace panel  = "`panel_name'"              in `row'
replace effect = "Proportion mediated (%)"   in `row'
local pm3 = prop_mediated3_std*100
if !missing(`pm3') replace pm_pct = string(`pm3',"%6.3f") in `row'

order panel effect beta_se ci95 p_value pm_pct
export delimited using "output/tables/mediation_analysis_results.csv", replace
restore

* ==============================================================================
* Clean up temporary scalars
* ==============================================================================

capture scalar drop a1_std b1_std c_total1_std c_direct1_std indirect1_std ///
    se_a1_std se_b1_std se_c_total1_std se_c_direct1_std se_indirect1_std ///
    p_a1_std p_b1_std p_c_total1_std p_c_direct1_std p_indirect1_std ///
    z_indirect1_std prop_mediated1_std z_c_total1_std z_c_direct1_std

capture scalar drop a2_std b2_std c_total2_std c_direct2_std indirect2_std ///
    se_a2_std se_b2_std se_c_total2_std se_c_direct2_std se_indirect2_std ///
    p_a2_std p_b2_std p_c_total2_std p_c_direct2_std p_indirect2_std ///
    z_indirect2_std prop_mediated2_std z_c_total2_std z_c_direct2_std

capture scalar drop a3_std b3_std c_total3_std c_direct3_std indirect3_std ///
    se_a3_std se_b3_std se_c_total3_std se_c_direct3_std se_indirect3_std ///
    p_a3_std p_b3_std p_c_total3_std p_c_direct3_std p_indirect3_std ///
    z_indirect3_std prop_mediated3_std z_c_total3_std z_c_direct3_std

di " "
di "=== All analyses completed. ==="
