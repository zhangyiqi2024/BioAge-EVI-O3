* =============================================================================
* Robust Imputation Scheme - Single Imputation with IPW Weights
* =============================================================================
set seed 20240101

version 17
clear all
set more off

* Load analysis-ready dataset
use "DATA/analysis_ready/analysis_dataset.dta", clear

* -----------------------------------------------------------------------------
* Step 1: Create trimmed variable set
* -----------------------------------------------------------------------------
di as txt "=== Creating trimmed variable set ==="

* Generate exposure lag-mean variables using loop
foreach prefix in EVI EVI_500m EVI_1KM O3 O3_500m O3_1KM {
    gen `prefix'_mean2y = (`prefix'_lag1 + `prefix'_lag2) / 2
    gen `prefix'_mean3y = (`prefix'_lag1 + `prefix'_lag2 + `prefix'_lag3) / 3
    gen `prefix'_mean4y = (`prefix'_lag1 + `prefix'_lag2 + `prefix'_lag3 + `prefix'_lag4) / 4
    gen `prefix'_mean5y = (`prefix'_lag1 + `prefix'_lag2 + `prefix'_lag3 + `prefix'_lag4 + `prefix'_lag5) / 5
}

* Generate categorical variables for O3_1KM_mean4y and EVI_1KM_mean4y
gen O3_1KM_mean4y_cat = cond(O3_1KM_mean4y < 85.85, 0, 1) if !missing(O3_1KM_mean4y)
gen EVI_1KM_mean4y_cat = cond(EVI_1KM_mean4y < 0.4, 0, 1) if !missing(EVI_1KM_mean4y)

* Set variable labels and value labels
label variable O3_1KM_mean4y_cat "O3_1KM_mean4y category"
label variable EVI_1KM_mean4y_cat "EVI_1KM_mean4y category"

label define o3_cat_label 0 "O3_1KM_mean4y<85.85" 1 "O3_1KM_mean4y>=85.85"
label define evi_cat_label 0 "EVI_1KM_mean4y<0.4" 1 "EVI_1KM_mean4y>=0.4"

label values O3_1KM_mean4y_cat o3_cat_label
label values EVI_1KM_mean4y_cat evi_cat_label

* Generate biological age acceleration binary variable using cond()
gen ba_acceleration_A = cond(BAA_PhenoAge <= 0, 0, 1) if !missing(BAA_PhenoAge)

* Set variable labels
label variable ba_acceleration_A "Biological age acceleration category"

* Set value labels
label define bioage_label 0 "BAA_PhenoAge<=0" 1 "BAA_PhenoAge>0"
label values ba_acceleration_A bioage_label

local wanted id city wave year ///
    hukou_type hh_income_total hh_consumption_total hh_consumption_percapita ///
    retirement pension_insurance medical_insurance smoking_current alcohol_current hobbies employed han_ethnicity ///
    hypertension dyslipidemia arthritis disability asthma sleep_quality hospital urban_rural ///
    bmi outpatient_cost_total cognitive_ability metabolic_activity mental_health bl_hbalc bl_tg bl_hdl bl_ldl ///
    resilience_A resilience_B TMP PER SNOW DEM lat WIND ///
    age sex marital_status adl_difficulty education_raw education self_rated_health life_satisfaction ///
    social_isolation1 social_isolation2 social_participation_score ///
    outpatient_visits inpatient_visits ///
    delta_BAA_PhenoAge BAA_PhenoAge BA_PhenoAge BAA_PhenoAge_base_2011 ///
    BA_PhenoAge_base_2011 ///
    pa_light pa_moderate pa_vigorous pa_light_moderate_weekly ///
    pa_moderate_vigorous_weekly pa_any_weekly social_participation cancer diabetes ///
    heart_disease stroke lung_disease memory_disease ///
    days_light days_moderate days_vigorous days_light_moderate days_moderate_vigorous ///
    days_any_pa ///
    traj15_NDVI traj15_NDVI_500m traj15_NDVI_1KM ///
    traj15_FVC traj15_FVC_500m traj15_FVC_1KM ///
    traj15_EVI traj15_EVI_500m traj15_EVI_1KM ///
    traj15_O3 traj15_O3_500m traj15_O3_1KM	 ///
    traj15_PM25 traj15_PM25_500m traj15_PM25_1KM ///
    traj15_PM10 traj15_PM10_500m traj15_PM10_1KM ///
    PM25_1KM PM10_1KM  ///
    PM25_500m PM10_500m ///
    PM25 PM10 ///
    PM25_1KM_lag1 PM10_1KM_lag1  ///
    PM25_500m_lag1 PM10_500m_lag1  ///
    PM25_lag1 PM10_lag1  ///
    EVI      EVI_500m      EVI_1KM  ///
    EVI_lag1 EVI_500m_lag1 EVI_1KM_lag1  ///
    EVI_lag2 EVI_500m_lag2 EVI_1KM_lag2  ///
    EVI_lag3 EVI_500m_lag3 EVI_1KM_lag3  ///
    EVI_lag4 EVI_500m_lag4 EVI_1KM_lag4  ///
    EVI_lag5 EVI_500m_lag5 EVI_1KM_lag5  ///
    O3       O3_500m       O3_1KM   ///
    O3_lag1  O3_500m_lag1  O3_1KM_lag1   ///
    O3_lag2  O3_500m_lag2  O3_1KM_lag2   ///
    O3_lag3  O3_500m_lag3  O3_1KM_lag3   ///
    O3_lag4  O3_500m_lag4  O3_1KM_lag4   ///
    O3_lag5  O3_500m_lag5  O3_1KM_lag5   ///
    EVI_mean2y EVI_mean3y EVI_mean4y EVI_mean5y   ///
    EVI_500m_mean2y EVI_500m_mean3y EVI_500m_mean4y EVI_500m_mean5y   ///
    EVI_1KM_mean2y EVI_1KM_mean3y EVI_1KM_mean4y EVI_1KM_mean5y   ///
    O3_mean2y O3_mean3y O3_mean4y O3_mean5y   ///
    O3_500m_mean2y O3_500m_mean3y O3_500m_mean4y O3_500m_mean5y   ///
    O3_1KM_mean2y O3_1KM_mean3y O3_1KM_mean4y O3_1KM_mean5y   ///
    O3_1KM_mean4y_cat EVI_1KM_mean4y_cat ba_acceleration_A  ///
    inversion_days1_lag5 inversion_days2_lag5 ///
    steep_slope_pct_25deg_1km_lag5 steep_slope_pct_25deg_500m_lag5 steep_slope_pct_25deg_lag5 ///
    steep_slope_pct_20deg_1km_lag5 steep_slope_pct_20deg_500m_lag5 steep_slope_pct_20deg_lag5 ///
    mean_slope_1km_lag5 mean_slope_500m_lag5 mean_slope_lag5 ///
    region region_class1 region_class2 hu_line    ///
    province_code city_code

local keep_list
foreach v of local wanted {
    capture confirm variable `v'
    if !_rc {
        local keep_list "`keep_list' `v'"
    }
    else {
        di as err "Warning: Variable `v' does not exist, skipped"
    }
}
keep `keep_list'

* Create individual identifier (if not present)
capture confirm variable id_num
if _rc {
    egen id_num = group(id)
    di as txt "Created individual identifier id_num"
}

save "DATA/analysis_ready/A_data_bio_age_PhenoAge.dta", replace

* Check for missing values in selected variables
local vars_to_check region_class1 region_class2 hu_line PM25_1KM PM10_1KM O3_1KM PM25_500m PM10_500m O3_500m PM25 PM10 O3 PM25_1KM_lag1 PM10_1KM_lag1 O3_1KM_lag1 PM25_500m_lag1 PM10_500m_lag1 O3_500m_lag1 PM25_lag1 PM10_lag1 O3_lag1

foreach var of local vars_to_check {
    capture confirm variable `var'
    if _rc == 0 {
        quietly count if missing(`var')
        if r(N) > 0 di "`var' has missing values: " r(N)
    }
}


* -----------------------------------------------------------------------------
* Step 2: Identify variables with missing values
* -----------------------------------------------------------------------------
di as txt "=== Identifying variables with missing values ==="

local all_vars ///
    hukou_type hh_income_total hh_consumption_total hh_consumption_percapita ///
    retirement pension_insurance smoking_current alcohol_current ///
    TMP PER WIND DEM lat disability ///
    han_ethnicity employed hypertension dyslipidemia diabetes cancer lung_disease heart_disease stroke memory_disease ///
    arthritis asthma pa_vigorous pa_moderate pa_light hobbies hospital ///
    medical_insurance social_participation pa_moderate_vigorous_weekly ///
    outpatient_visits inpatient_visits days_vigorous days_moderate days_light ///
    self_rated_health adl_difficulty ///
    cognitive_ability mental_health life_satisfaction metabolic_activity bmi outpatient_cost_total ///
    resilience_A resilience_B SNOW PM25_1KM PM10_1KM O3_1KM PM25_500m PM10_500m O3_500m  ///
    PM25 PM10 O3 PM25_1KM_lag1 PM10_1KM_lag1 O3_1KM_lag1  ///
    PM25_500m_lag1 PM10_500m_lag1 O3_500m_lag1 PM25_lag1 PM10_lag1 O3_lag1



local vars_to_impute
foreach var of local all_vars {
    capture confirm variable `var'
    if !_rc {
        quietly count if missing(`var')
        if r(N) > 0 {
            local vars_to_impute "`vars_to_impute' `var'"
            di as txt "Variable " as res "`var'" as txt " has " as res r(N) as txt " missing values"
        }
    }
}

di as txt "Total variables to impute: " as res `: word count `vars_to_impute''

if "`vars_to_impute'" == "" {
    di as txt "No variables need imputation, program ends"
    exit
}

* -----------------------------------------------------------------------------
* Step 3: Create backup of original data
* -----------------------------------------------------------------------------
foreach var of local vars_to_impute {
    capture confirm variable `var'
    if !_rc {
        generate `var'_orig = `var'
    }
}

* -----------------------------------------------------------------------------
* Step 4: Robust single imputation
* -----------------------------------------------------------------------------
di as txt "=== Starting robust single imputation ==="

foreach var of local vars_to_impute {
    capture confirm numeric variable `var'
    if !_rc {
        quietly count if missing(`var')
        local missing_count = r(N)
        if `missing_count' > 0 {
            di as txt "Imputing variable: " as res "`var'" as txt " (missing " as res `missing_count' as txt ")..."

            * Check variable type and distribution
            quietly summarize `var'
            local vmin = r(min)
            local vmax = r(max)
            local vmean = r(mean)
            local vsd = r(sd)
            local vN = r(N)

            * For binary variables (0/1)
            if `vmin' == 0 & `vmax' == 1 {
                di as txt "  - Type: Binary, proportion: " as res %6.4f `vmean'

                * Predict using logistic regression
                capture logit `var' urban_rural education age sex marital_status traj15_EVI_500m traj15_O3_500m BAA_PhenoAge
                if !_rc {
                    predict temp_pred
                    replace `var' = (temp_pred > 0.5) if missing(`var')
                    drop temp_pred
                    di as txt "  - Method: Logistic regression prediction"
                }
                else {
                    * If regression fails, use observed proportion for random imputation
                    local prob = `vmean'
                    replace `var' = (runiform() < `prob') if missing(`var')
                    di as txt "  - Method: Random imputation by observed proportion"
                }
            }
            * For count variables (non-negative integers)
            else if `vmin' >= 0 & `vmax' == round(`vmax') & `vmax' > 1 {
                di as txt "  - Type: Count, range: " as res `vmin' as txt "-" as res `vmax'

                * Predict using Poisson regression
                capture poisson `var' urban_rural education age sex marital_status traj15_EVI_500m traj15_O3_500m BAA_PhenoAge
                if !_rc {
                    predict temp_pred
                    replace `var' = round(temp_pred) if missing(`var')
                    replace `var' = max(`vmin', min(`var', `vmax')) if missing(`var'_orig)  // Ensure within valid range
                    drop temp_pred
                    di as txt "  - Method: Poisson regression prediction"
                }
                else {
                    * If regression fails, use observed mean for random imputation
                    local lambda = `vmean'
                    replace `var' = round(rpoisson(`lambda')) if missing(`var')
                    replace `var' = max(`vmin', min(`var', `vmax')) if missing(`var'_orig)
                    di as txt "  - Method: Random imputation from Poisson distribution"
                }
            }
            * For continuous variables
            else {
                di as txt "  - Type: Continuous, mean: " as res %8.2f `vmean' as txt ", SD: " as res %8.2f `vsd'

                * Predict using linear regression
                capture regress `var' urban_rural education age sex marital_status traj15_EVI_500m traj15_O3_500m BAA_PhenoAge
                if !_rc {
                    predict temp_pred
                    replace `var' = temp_pred if missing(`var')
                    drop temp_pred
                    di as txt "  - Method: Linear regression prediction"
                }
                else {
                    * If regression fails, use mean + random noise
                    replace `var' = `vmean' + rnormal(0, `vsd') if missing(`var')
                    di as txt "  - Method: Mean + random noise imputation"
                }
            }

            * Check post-imputation statistics
            quietly summarize `var'
            di as txt "  - After imputation: Mean = " as res %8.4f r(mean) as txt ", Missing = " as res r(N)-r(N_nonmiss)
        }
    }
}

* -----------------------------------------------------------------------------
* Step 5: Compute inverse probability weighting (IPW) for each core exposure variable (panel data version)
* -----------------------------------------------------------------------------
di as txt "=== Computing IPW weights for core exposure variables (panel data version) ==="

* Define core exposure variables (time-invariant)
local core_vars traj15_EVI_500m traj15_EVI_1KM traj15_O3_500m traj15_O3_1KM O3_1KM_mean4y_cat EVI_1KM_mean4y_cat

* Check whether core variables exist and are time-invariant
foreach var of local core_vars {
    capture confirm variable `var'
    if !_rc {
        di as txt "Checking time invariance of variable `var'..."

        * Check whether variable changes over time within each individual
        bysort id_num: egen var_sd = sd(`var')
        summarize var_sd
        if r(max) > 0.0001 {
            di as err "Warning: Variable `var' varies within individuals, violating panel data assumption"
        }
        else {
            di as txt "Variable `var' is time-invariant within individuals, consistent with panel data assumption"
        }
        drop var_sd

        di as txt "Computing IPW weights for `var'..."

        * Generate temporary variable name
        local ipw_var = "ipw_" + subinstr("`var'", "traj15_", "", 1)

        * Use baseline data (2011) to compute propensity scores
        preserve
        keep if year == 2011  // Use baseline data for weight computation

        * Check whether variable is binary
        quietly tab `var'
        if r(r) == 2 {
            * Use probit model to compute propensity scores (based on baseline covariates)
            tempvar phat
            capture probit `var' age sex education urban_rural marital_status
            if !_rc {
                predict double `phat', pr
                replace `phat' = 0.01 if `phat' < 0.01
                replace `phat' = 0.99 if `phat' > 0.99

                generate double `ipw_var' = 1/`phat' if `var' == 1 & !missing(`phat')
                replace `ipw_var' = 1/(1-`phat') if `var' == 0 & !missing(`phat')

                * Winsorize weights to avoid extreme values
                summarize `ipw_var', detail
                replace `ipw_var' = r(p1) if `ipw_var' < r(p1) & !missing(`ipw_var')
                replace `ipw_var' = r(p99) if `ipw_var' > r(p99) & !missing(`ipw_var')

                di as txt "IPW weight `ipw_var' statistics (baseline data):"
                di as txt "  Mean = " as res %6.4f r(mean)
                di as txt "  SD   = " as res %6.4f r(sd)
                di as txt "  Min  = " as res %6.4f r(min)
                di as txt "  Max  = " as res %6.4f r(max)

                * Save weights and id for merging
                keep id_num `ipw_var'
                save "temp_ipw_`var'.dta", replace
            }
            else {
                di as err "Propensity score model failed, generating equal weights"
                generate double `ipw_var' = 1
                keep id_num `ipw_var'
                save "temp_ipw_`var'.dta", replace
            }
        }
        else {
            di as err "Variable `var' is not binary, cannot compute IPW weights"
            generate double `ipw_var' = 1
            keep id_num `ipw_var'
            save "temp_ipw_`var'.dta", replace
        }
        restore

        * Merge weights back to full dataset
        merge m:1 id_num using "temp_ipw_`var'.dta", nogen
        erase "temp_ipw_`var'.dta"

        di as txt "Completed IPW weight computation and merge for variable `var'"
    }
    else {
        di as err "Variable `var' does not exist, skipping"
    }
}
* -----------------------------------------------------------------------------
* Step 5a: Check and compute missing single-exposure weights
* -----------------------------------------------------------------------------
di as txt "=== Checking and supplementing missing single-exposure weights ==="

* Check whether ipw_EVI_1KM_mean4y_cat has been computed
capture confirm variable ipw_EVI_1KM_mean4y_cat
if _rc {
    di as txt "ipw_EVI_1KM_mean4y_cat weight not found, computing now..."

    preserve
    keep if year == 2011  // Use baseline data for weight computation

    * Check whether variable is binary
    capture confirm variable EVI_1KM_mean4y_cat
    if !_rc {
        quietly tab EVI_1KM_mean4y_cat
        if r(r) == 2 {
            * Use probit model to compute propensity scores
            tempvar phat_evi_mean4y
            capture probit EVI_1KM_mean4y_cat age sex education urban_rural marital_status
            if !_rc {
                predict double `phat_evi_mean4y', pr
                replace `phat_evi_mean4y' = 0.01 if `phat_evi_mean4y' < 0.01
                replace `phat_evi_mean4y' = 0.99 if `phat_evi_mean4y' > 0.99

                generate double ipw_EVI_1KM_mean4y_cat = 1/`phat_evi_mean4y' if EVI_1KM_mean4y_cat == 1 & !missing(`phat_evi_mean4y')
                replace ipw_EVI_1KM_mean4y_cat = 1/(1-`phat_evi_mean4y') if EVI_1KM_mean4y_cat == 0 & !missing(`phat_evi_mean4y')

                * Winsorize weights
                summarize ipw_EVI_1KM_mean4y_cat, detail
                replace ipw_EVI_1KM_mean4y_cat = r(p1) if ipw_EVI_1KM_mean4y_cat < r(p1) & !missing(ipw_EVI_1KM_mean4y_cat)
                replace ipw_EVI_1KM_mean4y_cat = r(p99) if ipw_EVI_1KM_mean4y_cat > r(p99) & !missing(ipw_EVI_1KM_mean4y_cat)

                di as txt "IPW weight ipw_EVI_1KM_mean4y_cat statistics (baseline data):"
                di as txt "  Mean = " as res %6.4f r(mean)
                di as txt "  SD   = " as res %6.4f r(sd)
                di as txt "  Min  = " as res %6.4f r(min)
                di as txt "  Max  = " as res %6.4f r(max)

                * Save weights
                keep id_num ipw_EVI_1KM_mean4y_cat
                save "temp_ipw_EVI_1KM_mean4y_cat.dta", replace
            }
            else {
                di as err "Propensity score model failed, generating equal weights"
                generate double ipw_EVI_1KM_mean4y_cat = 1
                keep id_num ipw_EVI_1KM_mean4y_cat
                save "temp_ipw_EVI_1KM_mean4y_cat.dta", replace
            }
        }
        else {
            di as err "Variable EVI_1KM_mean4y_cat is not binary, cannot compute IPW weights"
            generate double ipw_EVI_1KM_mean4y_cat = 1
            keep id_num ipw_EVI_1KM_mean4y_cat
            save "temp_ipw_EVI_1KM_mean4y_cat.dta", replace
        }
    }
    else {
        di as err "Variable EVI_1KM_mean4y_cat does not exist"
    }
    restore

    * Merge weights
    capture confirm file "temp_ipw_EVI_1KM_mean4y_cat.dta"
    if !_rc {
        merge m:1 id_num using "temp_ipw_EVI_1KM_mean4y_cat.dta", nogen
        erase "temp_ipw_EVI_1KM_mean4y_cat.dta"
        di as txt "Completed IPW weight computation and merge for EVI_1KM_mean4y_cat"
    }
    else {
        di as err "temp_ipw_EVI_1KM_mean4y_cat.dta file not found, weight computation failed"
    }
}
else {
    di as txt "ipw_EVI_1KM_mean4y_cat weight already exists, skipping"
}

* Check whether ipw_O3_1KM_mean4y_cat has been computed (re-confirm)
capture confirm variable ipw_O3_1KM_mean4y_cat
if _rc {
    di as txt "ipw_O3_1KM_mean4y_cat weight not found, computing now..."

    preserve
    keep if year == 2011

    capture confirm variable O3_1KM_mean4y_cat
    if !_rc {
        quietly tab O3_1KM_mean4y_cat
        if r(r) == 2 {
            tempvar phat_o3_mean4y
            capture probit O3_1KM_mean4y_cat age sex education urban_rural marital_status
            if !_rc {
                predict double `phat_o3_mean4y', pr
                replace `phat_o3_mean4y' = 0.01 if `phat_o3_mean4y' < 0.01
                replace `phat_o3_mean4y' = 0.99 if `phat_o3_mean4y' > 0.99

                generate double ipw_O3_1KM_mean4y_cat = 1/`phat_o3_mean4y' if O3_1KM_mean4y_cat == 1 & !missing(`phat_o3_mean4y')
                replace ipw_O3_1KM_mean4y_cat = 1/(1-`phat_o3_mean4y') if O3_1KM_mean4y_cat == 0 & !missing(`phat_o3_mean4y')

                summarize ipw_O3_1KM_mean4y_cat, detail
                replace ipw_O3_1KM_mean4y_cat = r(p1) if ipw_O3_1KM_mean4y_cat < r(p1) & !missing(ipw_O3_1KM_mean4y_cat)
                replace ipw_O3_1KM_mean4y_cat = r(p99) if ipw_O3_1KM_mean4y_cat > r(p99) & !missing(ipw_O3_1KM_mean4y_cat)

                di as txt "IPW weight ipw_O3_1KM_mean4y_cat statistics (baseline data):"
                di as txt "  Mean = " as res %6.4f r(mean)
                di as txt "  SD   = " as res %6.4f r(sd)

                keep id_num ipw_O3_1KM_mean4y_cat
                save "temp_ipw_O3_1KM_mean4y_cat.dta", replace
            }
            else {
                generate double ipw_O3_1KM_mean4y_cat = 1
                keep id_num ipw_O3_1KM_mean4y_cat
                save "temp_ipw_O3_1KM_mean4y_cat.dta", replace
            }
        }
        else {
            generate double ipw_O3_1KM_mean4y_cat = 1
            keep id_num ipw_O3_1KM_mean4y_cat
            save "temp_ipw_O3_1KM_mean4y_cat.dta", replace
        }
    }
    restore

    capture confirm file "temp_ipw_O3_1KM_mean4y_cat.dta"
    if !_rc {
        merge m:1 id_num using "temp_ipw_O3_1KM_mean4y_cat.dta", nogen
        erase "temp_ipw_O3_1KM_mean4y_cat.dta"
        di as txt "Completed IPW weight computation and merge for O3_1KM_mean4y_cat"
    }
}
else {
    di as txt "ipw_O3_1KM_mean4y_cat weight already exists"
}

di as txt "Single-exposure weight check and supplement complete"
* -----------------------------------------------------------------------------
* Step 5b: Compute inverse probability weights based on joint exposure classification (following literature approach)
* -----------------------------------------------------------------------------
di as txt "=== Computing IPW based on joint exposure classification (4 exposure groups) ==="

* Check whether exposure variables are binary
foreach var in traj15_EVI_1KM traj15_EVI_500m traj15_O3_1KM traj15_O3_500m O3_1KM_mean4y_cat EVI_1KM_mean4y_cat {
    capture confirm variable `var'
    if !_rc {
        quietly tab `var'
        di as txt "Variable `var': " as res r(r) as txt " categories"
        if r(r) != 2 {
            di as err "Warning: Variable `var' is not binary; joint exposure weight computation may not apply"
        }
    }
}

* 1. Compute joint exposure weights at 1KM scale
di as txt "Computing joint exposure weights at 1KM scale..."

preserve
keep if year == 2011  // Use baseline data for weight computation

* Create joint exposure classification (4 categories)
capture confirm variable traj15_EVI_1KM
capture confirm variable traj15_O3_1KM
if !_rc {
    * Generate joint exposure classification (4 categories), following literature approach
    generate joint_exp_1KM = .
    replace joint_exp_1KM = 1 if traj15_EVI_1KM == 0 & traj15_O3_1KM == 0  // Low EVI, Low O3
    replace joint_exp_1KM = 2 if traj15_EVI_1KM == 0 & traj15_O3_1KM == 1  // Low EVI, High O3
    replace joint_exp_1KM = 3 if traj15_EVI_1KM == 1 & traj15_O3_1KM == 0  // High EVI, Low O3
    replace joint_exp_1KM = 4 if traj15_EVI_1KM == 1 & traj15_O3_1KM == 1  // High EVI, High O3

    label define joint_exp_1KM 1 "Low EVI, Low O3" 2 "Low EVI, High O3" 3 "High EVI, Low O3" 4 "High EVI, High O3"
    label values joint_exp_1KM joint_exp_1KM

    * Check joint exposure classification distribution
    di as txt "1KM-scale joint exposure classification distribution:"
    tab joint_exp_1KM, mi

    * Use multinomial logit model to compute joint propensity scores
    tempvar phat_joint1
    capture mlogit joint_exp_1KM age sex education urban_rural marital_status
    if !_rc {
        predict double `phat_joint1'

        * Compute inverse probability weight: 1/predicted probability
        generate double ipw_joint_1KM = 1/`phat_joint1'

        * Winsorize weights (1%-99%)
        summarize ipw_joint_1KM, detail
        replace ipw_joint_1KM = r(p1) if ipw_joint_1KM < r(p1) & !missing(ipw_joint_1KM)
        replace ipw_joint_1KM = r(p99) if ipw_joint_1KM > r(p99) & !missing(ipw_joint_1KM)

        di as txt "Joint weight ipw_joint_1KM statistics:"
        di as txt "  Mean = " as res %6.4f r(mean)
        di as txt "  SD   = " as res %6.4f r(sd)
        di as txt "  Min  = " as res %6.4f r(min)
        di as txt "  Max  = " as res %6.4f r(max)

        * Check weight distribution by exposure category
        di as txt "Weight distribution by joint exposure category:"
        tabstat ipw_joint_1KM, by(joint_exp_1KM) stat(mean sd min max n) nototal
    }
    else {
        di as err "Multinomial logit model failed, using product method for joint weights"
        generate double ipw_joint_1KM = ipw_EVI_1KM * ipw_O3_1KM
    }

    * Save weights and joint exposure classification
    keep id_num ipw_joint_1KM joint_exp_1KM
    save "temp_joint_ipw_1KM.dta", replace
}
restore

* Merge 1KM joint weights
capture confirm file "temp_joint_ipw_1KM.dta"
if !_rc {
    merge m:1 id_num using "temp_joint_ipw_1KM.dta", nogen
    erase "temp_joint_ipw_1KM.dta"
    di as txt "1KM-scale joint exposure weights merged successfully"
}

*-----------------------------------------------------------------------

* 2. Compute joint exposure weights at 500m scale
di as txt "Computing joint exposure weights at 500m scale..."

preserve
keep if year == 2011  // Use baseline data for weight computation

* Create joint exposure classification (4 categories)
capture confirm variable traj15_EVI_500m
capture confirm variable traj15_O3_500m
if !_rc {
    * Generate joint exposure classification (4 categories), following literature approach
    generate joint_exp_500m = .
    replace joint_exp_500m = 1 if traj15_EVI_500m == 0 & traj15_O3_500m == 0  // Low EVI, Low O3
    replace joint_exp_500m = 2 if traj15_EVI_500m == 0 & traj15_O3_500m == 1  // Low EVI, High O3
    replace joint_exp_500m = 3 if traj15_EVI_500m == 1 & traj15_O3_500m == 0  // High EVI, Low O3
    replace joint_exp_500m = 4 if traj15_EVI_500m == 1 & traj15_O3_500m == 1  // High EVI, High O3

    label define joint_exp_500m 1 "Low EVI, Low O3" 2 "Low EVI, High O3" 3 "High EVI, Low O3" 4 "High EVI, High O3"
    label values joint_exp_500m joint_exp_500m

    * Check joint exposure classification distribution
    di as txt "500m-scale joint exposure classification distribution:"
    tab joint_exp_500m, mi

    * Use multinomial logit model to compute joint propensity scores
    tempvar phat_joint2
    capture mlogit joint_exp_500m age sex education urban_rural marital_status
    if !_rc {
        predict double `phat_joint2'

        * Compute inverse probability weight: 1/predicted probability
        generate double ipw_joint_500m = 1/`phat_joint2'

        * Winsorize weights (1%-99%)
        summarize ipw_joint_500m, detail
        replace ipw_joint_500m = r(p1) if ipw_joint_500m < r(p1) & !missing(ipw_joint_500m)
        replace ipw_joint_500m = r(p99) if ipw_joint_500m > r(p99) & !missing(ipw_joint_500m)

        di as txt "Joint weight ipw_joint_500m statistics:"
        di as txt "  Mean = " as res %6.4f r(mean)
        di as txt "  SD   = " as res %6.4f r(sd)
        di as txt "  Min  = " as res %6.4f r(min)
        di as txt "  Max  = " as res %6.4f r(max)

        * Check weight distribution by exposure category
        di as txt "Weight distribution by joint exposure category:"
        tabstat ipw_joint_500m, by(joint_exp_500m) stat(mean sd min max n) nototal
    }
    else {
        di as err "Multinomial logit model failed, using product method for joint weights"
        generate double ipw_joint_500m = ipw_EVI_500m * ipw_O3_500m
    }

    * Save weights and joint exposure classification
    keep id_num ipw_joint_500m joint_exp_500m
    save "temp_joint_ipw_500m.dta", replace
}
restore

* Merge 500m joint weights
capture confirm file "temp_joint_ipw_500m.dta"
if !_rc {
    merge m:1 id_num using "temp_joint_ipw_500m.dta", nogen
    erase "temp_joint_ipw_500m.dta"
    di as txt "500m-scale joint exposure weights merged successfully"
}

di as txt "Joint exposure classification weight computation complete"

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

* 3. Compute joint exposure weights for EVI_1KM_mean4y_cat x O3_1KM_mean4y_cat at 1KM scale
* Note: These exposure variables are time-varying, but we use baseline (2011) exposure status for weights
di as txt "Computing 1KM-scale mean joint exposure weights (based on baseline exposure status)..."

* Check whether necessary single weights exist
capture confirm variable ipw_EVI_1KM_mean4y_cat
if _rc {
    di as err "Error: ipw_EVI_1KM_mean4y_cat weight does not exist"
    exit
}

capture confirm variable ipw_O3_1KM_mean4y_cat
if _rc {
    di as err "Error: ipw_O3_1KM_mean4y_cat weight does not exist"
    exit
}

preserve
* Key modification: Use only baseline (2011) data to compute joint exposure classification and weights
keep if year == 2011

* Check whether exposure variables exist
capture confirm variable EVI_1KM_mean4y_cat
capture confirm variable O3_1KM_mean4y_cat
if !_rc {
    * Drop potentially existing temporary variables
    capture drop joint_exp_1KM_mean4y_temp ipw_joint_1KM_mean4y_temp

    * Generate baseline joint exposure classification (4 categories)
    generate joint_exp_1KM_mean4y_temp = .
    replace joint_exp_1KM_mean4y_temp = 1 if EVI_1KM_mean4y_cat == 0 & O3_1KM_mean4y_cat == 0  // Low EVI, Low O3
    replace joint_exp_1KM_mean4y_temp = 2 if EVI_1KM_mean4y_cat == 0 & O3_1KM_mean4y_cat == 1  // Low EVI, High O3
    replace joint_exp_1KM_mean4y_temp = 3 if EVI_1KM_mean4y_cat == 1 & O3_1KM_mean4y_cat == 0  // High EVI, Low O3
    replace joint_exp_1KM_mean4y_temp = 4 if EVI_1KM_mean4y_cat == 1 & O3_1KM_mean4y_cat == 1  // High EVI, High O3

    label define joint_exp_1KM_mean4y_temp 1 "Low EVI, Low O3" 2 "Low EVI, High O3" 3 "High EVI, Low O3" 4 "High EVI, High O3"
    label values joint_exp_1KM_mean4y_temp joint_exp_1KM_mean4y_temp

    * Check baseline joint exposure classification distribution
    di as txt "Baseline (2011) 1KM-scale mean joint exposure classification distribution:"
    tab joint_exp_1KM_mean4y_temp, mi

    * Use multinomial logit model to compute joint propensity scores (based on baseline covariates)
    tempvar phat_joint3
    capture mlogit joint_exp_1KM_mean4y_temp age sex education urban_rural marital_status
    if !_rc {
        predict double `phat_joint3'

        * Compute inverse probability weight: 1/predicted probability
        generate double ipw_joint_1KM_mean4y_temp = 1/`phat_joint3'

        * Winsorize weights (1%-99%)
        summarize ipw_joint_1KM_mean4y_temp, detail
        replace ipw_joint_1KM_mean4y_temp = r(p1) if ipw_joint_1KM_mean4y_temp < r(p1) & !missing(ipw_joint_1KM_mean4y_temp)
        replace ipw_joint_1KM_mean4y_temp = r(p99) if ipw_joint_1KM_mean4y_temp > r(p99) & !missing(ipw_joint_1KM_mean4y_temp)

        di as txt "Baseline joint weight ipw_joint_1KM_mean4y statistics:"
        di as txt "  Mean = " as res %6.4f r(mean)
        di as txt "  SD   = " as res %6.4f r(sd)
        di as txt "  Min  = " as res %6.4f r(min)
        di as txt "  Max  = " as res %6.4f r(max)

        * Check weight distribution by exposure category
        di as txt "Baseline weight distribution by joint exposure category:"
        tabstat ipw_joint_1KM_mean4y_temp, by(joint_exp_1KM_mean4y_temp) stat(mean sd min max n) nototal
    }
    else {
        di as err "Multinomial logit model failed, using product method for joint weights"
        * Use product of baseline single weights
        generate double ipw_joint_1KM_mean4y_temp = ipw_EVI_1KM_mean4y_cat * ipw_O3_1KM_mean4y_cat

        * Winsorize product weights
        summarize ipw_joint_1KM_mean4y_temp, detail
        replace ipw_joint_1KM_mean4y_temp = r(p1) if ipw_joint_1KM_mean4y_temp < r(p1) & !missing(ipw_joint_1KM_mean4y_temp)
        replace ipw_joint_1KM_mean4y_temp = r(p99) if ipw_joint_1KM_mean4y_temp > r(p99) & !missing(ipw_joint_1KM_mean4y_temp)

        di as txt "Product-method baseline joint weight ipw_joint_1KM_mean4y statistics:"
        di as txt "  Mean = " as res %6.4f r(mean)
        di as txt "  SD   = " as res %6.4f r(sd)
        di as txt "  Min  = " as res %6.4f r(min)
        di as txt "  Max  = " as res %6.4f r(max)
    }

    * Rename to final variable names
    rename joint_exp_1KM_mean4y_temp joint_exp_1KM_mean4y_base
    rename ipw_joint_1KM_mean4y_temp ipw_joint_1KM_mean4y_base

    * Save baseline weights and joint exposure classification
    keep id_num ipw_joint_1KM_mean4y_base joint_exp_1KM_mean4y_base
    save "temp_joint_ipw_1KM_mean4y_base.dta", replace

    di as txt "Baseline joint exposure weights saved"
}
else {
    di as err "EVI_1KM_mean4y_cat or O3_1KM_mean4y_cat variables do not exist"
}
restore

* Merge baseline joint weights to full dataset
capture confirm file "temp_joint_ipw_1KM_mean4y_base.dta"
if !_rc {
    * Drop potentially existing variables first
    capture drop joint_exp_1KM_mean4y_base ipw_joint_1KM_mean4y_base

    merge m:1 id_num using "temp_joint_ipw_1KM_mean4y_base.dta", nogen
    erase "temp_joint_ipw_1KM_mean4y_base.dta"

    * Create current exposure status joint classification for all observations (for descriptive analysis)
    generate joint_exp_1KM_mean4y_current = .
    replace joint_exp_1KM_mean4y_current = 1 if EVI_1KM_mean4y_cat == 0 & O3_1KM_mean4y_cat == 0
    replace joint_exp_1KM_mean4y_current = 2 if EVI_1KM_mean4y_cat == 0 & O3_1KM_mean4y_cat == 1
    replace joint_exp_1KM_mean4y_current = 3 if EVI_1KM_mean4y_cat == 1 & O3_1KM_mean4y_cat == 0
    replace joint_exp_1KM_mean4y_current = 4 if EVI_1KM_mean4y_cat == 1 & O3_1KM_mean4y_cat == 1

    label define joint_exp_1KM_mean4y_current 1 "Low EVI, Low O3" 2 "Low EVI, High O3" 3 "High EVI, Low O3" 4 "High EVI, High O3"
    label values joint_exp_1KM_mean4y_current joint_exp_1KM_mean4y_current

    di as txt "1KM-scale mean joint exposure weights merged successfully"
    di as txt "Note: Using baseline (2011) exposure status weights, applicable to entire follow-up period"
    di as txt "Current exposure status distribution:"
    tab joint_exp_1KM_mean4y_current year, mi
}
else {
    di as err "temp_joint_ipw_1KM_mean4y_base.dta file not found; joint weight computation failed"
}

* -----------------------------------------------------------------------------
* Step 6: Create analysis dataset
* -----------------------------------------------------------------------------

* Drop backup variables
capture drop *_orig

* Mark dataset type
generate dataset_type = "Robust single imputation, panel data"
generate imputation_date = "`c(current_date)'"
generate imputation_method = "Regression prediction based on age and sex"

* Ensure baseline BAA_PhenoAge variable exists
capture confirm variable BAA_PhenoAge_base_2011
if _rc {
    di as err "Warning: BAA_PhenoAge_base_2011 does not exist, creating it"
    bysort id_num: egen BAA_PhenoAge_base_2011 = min(BAA_PhenoAge) if year == 2011
    bysort id_num: egen temp_base = mean(BAA_PhenoAge_base_2011)
    replace BAA_PhenoAge_base_2011 = temp_base
    drop temp_base
}

* Save full dataset with all weight variables
save "DATA/analysis_ready/analysis_dataset.dta", replace
di as txt "Full dataset with weights saved"
