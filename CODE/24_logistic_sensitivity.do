* -----------------------------------------------------------------------------
* Mixed-effects Logistic Model (Random-effects Logistic Model)
* Method 1: Logistic Regression Model - Two Spatial Scales
* DV: ba_accelerated (based on BAA_PhenoAge quintiles: Q4-Q5=1, Q1-Q2=0, Q3 excluded)
* Exposure variables: EVI and O3 averaged over 4-year lag
*   - EVI: Enhanced Vegetation Index, raw scale 0-1, OR reported per 0.1 unit increase
*   - O3: Ozone concentration, raw scale 0-100, OR reported per 10 unit increase
* -----------------------------------------------------------------------------

set seed 20240101

di as txt "=== Method 1: Logistic Regression Model - Effect Size Comparison ==="

* Reload data
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

* -------------------------------------------------------------------------
* 1. Generate binary DV ba_accelerated from BAA_PhenoAge quintiles
* Steps:
* 1) xtile generates 5-quantile variable: baa_q = 1,2,3,4,5
* 2) baa_q >= 4 -> ba_accelerated = 1 (accelerated: Q4-Q5)
*    baa_q <= 2 -> ba_accelerated = 0 (not accelerated: Q1-Q2)
*    baa_q == 3 -> set to missing, excluded from regression
* -------------------------------------------------------------------------
capture drop baa_q
xtile baa_q = BAA_PhenoAge if !missing(BAA_PhenoAge), nq(5)

di as txt "=== BAA_PhenoAge Quintile Distribution Check ==="
tab baa_q

capture drop ba_accelerated
generate ba_accelerated = .
replace ba_accelerated = 1 if inlist(baa_q, 4, 5)
replace ba_accelerated = 0 if inlist(baa_q, 1, 2)
* baa_q == 3 remains missing, automatically excluded in logit

label define BAaccel 0 "Not accelerated (Q1-Q2)" 1 "Accelerated (Q4-Q5)"
label values ba_accelerated BAaccel

di as txt "=== ba_accelerated Variable Distribution Check (Q1-Q2 vs Q4-Q5 only) ==="
tab ba_accelerated

* -------------------------------------------------------------------------
* 2. Construct exposure variables scaled per 0.1 EVI and per 10 O3
*   EVI_..._p01: raw EVI / 0.1 -> coefficient corresponds to effect per 0.1 unit EVI
*   O3_..._p10: raw O3 / 10   -> coefficient corresponds to effect per 10 unit O3
* -------------------------------------------------------------------------
capture drop EVI_1KM_mean4y_p01 O3_1KM_mean4y_p10
capture drop EVI_500m_mean4y_p01 O3_500m_mean4y_p10

gen EVI_1KM_mean4y_p01   = EVI_1KM_mean4y / 0.1
label var EVI_1KM_mean4y_p01 "EVI_1KM_mean4y (per 0.1)"

gen O3_1KM_mean4y_p10    = O3_1KM_mean4y  / 10
label var O3_1KM_mean4y_p10 "O3_1KM_mean4y (per 10)"

gen EVI_500m_mean4y_p01  = EVI_500m_mean4y / 0.1
label var EVI_500m_mean4y_p01 "EVI_500m_mean4y (per 0.1)"

gen O3_500m_mean4y_p10   = O3_500m_mean4y  / 10
label var O3_500m_mean4y_p10 "O3_500m_mean4y (per 10)"

* Ensure estout package is installed
capture which esttab
if _rc != 0 {
    ssc install estout, replace
}

di "=== Climate Variable PCA ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* -----------------------------------------------------------------------------
* 3. 1KM-scale Logistic Regression Models
* DV: ba_accelerated
* Exposure: scaled EVI_1KM_mean4y_p01 (per 0.1) and O3_1KM_mean4y_p10 (per 10)
* -----------------------------------------------------------------------------

* 1KM-scale Logistic Regression - Joint model
di as txt "--- 1KM-scale Logistic Regression (EVI & O3 Joint) ---"
capture noisily logit ba_accelerated ///
     c.EVI_1KM_mean4y_p01 c.O3_1KM_mean4y_p10 ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year, ///
     vce(robust)

if _rc == 0 {
    estimates store logit_1km_joint
    di as txt "1KM joint model estimated successfully"
}
else {
    di as error "1KM joint model estimation failed, error code: " _rc
}

* 1KM-scale Logistic Regression - EVI-only model
di as txt "--- 1KM-scale Logistic Regression (EVI only) ---"
capture noisily logit ba_accelerated ///
     c.EVI_1KM_mean4y_p01 ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year, ///
     vce(robust)

if _rc == 0 {
    estimates store logit_1km_evi
    di as txt "1KM EVI-only model estimated successfully"
}
else {
    di as error "1KM EVI-only model estimation failed, error code: " _rc
}

* 1KM-scale Logistic Regression - O3-only model
di as txt "--- 1KM-scale Logistic Regression (O3 only) ---"
capture noisily logit ba_accelerated ///
     c.O3_1KM_mean4y_p10 ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year, ///
     vce(robust)

if _rc == 0 {
    estimates store logit_1km_o3
    di as txt "1KM O3-only model estimated successfully"
}
else {
    di as error "1KM O3-only model estimation failed, error code: " _rc
}

* -----------------------------------------------------------------------------
* 4. 500m-scale Logistic Regression Models
* Exposure: EVI_500m_mean4y_p01 (per 0.1) and O3_500m_mean4y_p10 (per 10)
* -----------------------------------------------------------------------------
di as txt "--- 500m-scale Logistic Regression (EVI & O3 Joint) ---"

* First check if 500m-scale variables exist
capture describe EVI_500m_mean4y O3_500m_mean4y
if _rc == 0 {
    di as txt "500m-scale variables exist, beginning model estimation"

    capture noisily logit ba_accelerated ///
         c.EVI_500m_mean4y_p01 c.O3_500m_mean4y_p10 ///
         c.age i.sex i.education ///
         c.ln_hh_consumption_total i.retirement i.pension_insurance ///
         i.smoking_current i.alcohol_current ///
         climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
         c.bmi i.hypertension i.dyslipidemia ///
         i.year, ///
         vce(robust)

    if _rc == 0 {
        estimates store logit_500m_joint
        di as txt "500m joint model estimated successfully"

        * Optional: display key results for 500m joint model (OR already per 0.1 EVI / per 10 O3)
        matrix b = e(b)
        matrix V = e(V)
        local EVI_coef = b[1,1]
        local EVI_se   = sqrt(V[1,1])
        local EVI_or   = exp(`EVI_coef')

        local O3_coef  = b[1,2]
        local O3_se    = sqrt(V[2,2])
        local O3_or    = exp(`O3_coef')

        di as txt "500m joint model key results (EVI per 0.1, O3 per 10):"
        di as txt "EVI_500m (per 0.1): OR = " as res %8.4f `EVI_or'
        di as txt "O3_500m  (per 10):  OR = " as res %8.4f `O3_or'
    }
    else {
        di as error "500m joint model estimation failed, error code: " _rc
    }
}
else {
    di as error "500m-scale variables not found, please check variable names"
    di as txt "Variable list in current dataset (containing EVI and O3):"
    ds *EVI* *O3*
}

* 500m-scale Logistic Regression - EVI-only model
di as txt "--- 500m-scale Logistic Regression (EVI only) ---"
capture noisily logit ba_accelerated ///
     c.EVI_500m_mean4y_p01 ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year, ///
     vce(robust)

if _rc == 0 {
    estimates store logit_500m_evi
    di as txt "500m EVI-only model estimated successfully"
}
else {
    di as error "500m EVI-only model estimation failed, error code: " _rc
}

* 500m-scale Logistic Regression - O3-only model
di as txt "--- 500m-scale Logistic Regression (O3 only) ---"
capture noisily logit ba_accelerated ///
     c.O3_500m_mean4y_p10 ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year, ///
     vce(robust)

if _rc == 0 {
    estimates store logit_500m_o3
    di as txt "500m O3-only model estimated successfully"
}
else {
    di as error "500m O3-only model estimation failed, error code: " _rc
}

* -----------------------------------------------------------------------------
* 5. Output summary table of all model results
*   - eform: output OR
*   - ci: output 95% CI (CI for OR)
*   - SE not displayed
*   - Note: EVI per 0.1 unit, O3 per 10 units
* -----------------------------------------------------------------------------
di as txt _newline "=== All Model Results Summary ==="

* Check which models were successfully estimated
di as txt "Stored models:"
estimates dir

* Use keep option to show only main exposure variables (already scaled per unit)
di as txt _newline "=== Main Variable Results Summary (All Models) ==="
capture noisily esttab logit_1km_joint logit_1km_evi logit_1km_o3 ///
                     logit_500m_joint logit_500m_evi logit_500m_o3, ///
    keep(EVI_1KM_mean4y_p01 O3_1KM_mean4y_p10 ///
         EVI_500m_mean4y_p01 O3_500m_mean4y_p10) ///
    eform b(3) ci(3) ///
    varwidth(30) ///
    title("Main Variable Results Summary - OR(95% CI); EVI per 0.1, O3 per 10") ///
    mtitle("1KM Joint" "1KM EVI-only" "1KM O3-only" "500m Joint" "500m EVI-only" "500m O3-only") ///
    compress
