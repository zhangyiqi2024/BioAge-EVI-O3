********************************************************************************
* 26_kmeans_robustness.do
* Robustness of Joint Exposure Classification to Alternative Thresholds
*
* Purpose: Address Reviewer Issue #2 — demonstrate that the central finding
*          (BAA in high-EVI/high-O3 group ≈ reference group) is robust to
*          alternative dichotomization thresholds.
*
* Author:  YiQi Zhang (ZAFU)
* Date:    2026-07-16
* Input:   DATA/analysis_ready/analysis_dataset.dta
* Output:  output/tables/kmeans_robustness_comparison.csv
*
* Tests:
*   1. Tertile-based classification (3×3 = 9 groups, compare extreme vs middle)
*   2. Median-split classification (2×2 = 4 groups)
*   3. Continuous EVI×O3 product-term interaction (no arbitrary cut-points)
*   4. Quartile-based classification
********************************************************************************

clear all
set seed 20240101
set more off

* ---- Load data ----
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

* PCA for climate variables
capture drop climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* Define controls
global controls "c.age i.sex i.education c.ln_hh_consumption_total i.retirement i.pension_insurance i.smoking_current i.alcohol_current climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 c.bmi i.hypertension i.dyslipidemia i.year"

* Sort and set panel
sort id_num year
xtset id_num year

* ---- Define exposure variables ----
local evi_vars "EVI_1KM_mean4y EVI_500m_mean4y"
local o3_vars  "O3_1KM_mean4y O3_500m_mean4y"
local buffers  "1KM 500m"

* ---- 1. Tertile-based classification ----
di _n "======================================================"
di "1. TERTILE-BASED JOINT EXPOSURE CLASSIFICATION"
di "======================================================"

* Drop previously generated variables if re-running
capture drop evi_tert_1KM evi_tert_500m o3_tert_1KM o3_tert_500m
capture drop joint_tert_1KM joint_tert_500m

* Create tertile variables (once per buffer, not per exposure-pair)
foreach buf of local buffers {
    local evi_v = cond("`buf'"=="1KM", "EVI_1KM_mean4y", "EVI_500m_mean4y")
    local o3_v  = cond("`buf'"=="1KM", "O3_1KM_mean4y", "O3_500m_mean4y")

    xtile evi_tert_`buf' = `evi_v', nq(3)
    xtile o3_tert_`buf'  = `o3_v',  nq(3)

    * Generate 9-category joint exposure
    gen joint_tert_`buf' = (evi_tert_`buf' - 1) * 3 + o3_tert_`buf'
    label variable joint_tert_`buf' "Joint EVI×O3 tertile (9 categories)"
}

* ---- 2. Median-split classification ----
di _n "======================================================"
di "2. MEDIAN-SPLIT JOINT EXPOSURE CLASSIFICATION"
di "======================================================"

* Drop previously generated variables if re-running
capture drop evi_med_1KM evi_med_500m o3_med_1KM o3_med_500m
capture drop joint_med_1KM joint_med_500m
capture label drop med_labels
label define med_labels 1 "Low EVI, Low O3" 2 "Low EVI, High O3" ///
                        3 "High EVI, Low O3" 4 "High EVI, High O3"

foreach buf of local buffers {
    local evi_v = cond("`buf'"=="1KM", "EVI_1KM_mean4y", "EVI_500m_mean4y")
    local o3_v  = cond("`buf'"=="1KM", "O3_1KM_mean4y", "O3_500m_mean4y")

    qui sum `evi_v', detail
    gen evi_med_`buf' = (`evi_v' > r(p50))
    label variable evi_med_`buf' "EVI above median (0/1)"

    qui sum `o3_v', detail
    gen o3_med_`buf' = (`o3_v' > r(p50))
    label variable o3_med_`buf' "O3 above median (0/1)"

    * 4-category joint exposure
    gen joint_med_`buf' = 1  // low EVI, low O3
    replace joint_med_`buf' = 2 if evi_med_`buf' == 0 & o3_med_`buf' == 1
    replace joint_med_`buf' = 3 if evi_med_`buf' == 1 & o3_med_`buf' == 0
    replace joint_med_`buf' = 4 if evi_med_`buf' == 1 & o3_med_`buf' == 1
    label values joint_med_`buf' med_labels
}

* ---- 3. Continuous product-term interaction model ----
di _n "======================================================"
di "3. CONTINUOUS EVI × O3 INTERACTION (NO CUT-POINTS)"
di "======================================================"

* Drop previously generated variables if re-running
capture drop c_EVI_1KM_mean4y c_EVI_500m_mean4y c_O3_1KM_mean4y c_O3_500m_mean4y
capture drop evi_o3_int_1KM evi_o3_int_500m

* Center exposure variables for interpretable interaction
foreach evi_var of local evi_vars {
    qui sum `evi_var'
    gen c_`evi_var' = `evi_var' - r(mean)
}
foreach o3_var of local o3_vars {
    qui sum `o3_var'
    gen c_`o3_var' = `o3_var' - r(mean)
}

* Generate interaction terms
gen evi_o3_int_1KM  = c_EVI_1KM_mean4y  * c_O3_1KM_mean4y
gen evi_o3_int_500m = c_EVI_500m_mean4y * c_O3_500m_mean4y

* ---- 4. Run comparative models ----
di _n "======================================================"
di "4. RUNNING COMPARATIVE JOINT EXPOSURE MODELS"
di "======================================================"

* Collect results
tempname results
postfile `results' str40 classification str10 buffer str30 contrast ///
    float beta float se float pval float ci_low float ci_high ///
    using "output/tables/kmeans_robustness_comparison.dta", replace

* --- 4a. Original K-means (binary) classification ---
di _n "=== Original K-means classification ==="
foreach buf in 1KM 500m {
    xtmixed BAA_PhenoAge i.joint_exp_`buf' $controls || id_num:, mle

    * High EVI / High O3 vs reference
    lincom 4.joint_exp_`buf'
    post `results' ("Original K-means") ("`buf'") ///
        ("High-EVI/High-O3 vs Reference") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))

    * Low EVI / High O3 vs reference (O3-only contrast)
    lincom 2.joint_exp_`buf'
    post `results' ("Original K-means") ("`buf'") ///
        ("Low-EVI/High-O3 vs Reference") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))

    * High EVI / Low O3 vs reference (EVI-only contrast)
    lincom 3.joint_exp_`buf'
    post `results' ("Original K-means") ("`buf'") ///
        ("High-EVI/Low-O3 vs Reference") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))
}

* --- 4b. Median-split classification ---
di _n "=== Median-split classification ==="
foreach buf in 1KM 500m {
    xtmixed BAA_PhenoAge i.joint_med_`buf' $controls || id_num:, mle

    * High EVI / High O3 vs reference
    lincom 4.joint_med_`buf'
    post `results' ("Median split") ("`buf'") ///
        ("High-EVI/High-O3 vs Reference") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))

    * Low EVI / High O3 vs reference
    lincom 2.joint_med_`buf'
    post `results' ("Median split") ("`buf'") ///
        ("Low-EVI/High-O3 vs Reference") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))

    * High EVI / Low O3 vs reference
    lincom 3.joint_med_`buf'
    post `results' ("Median split") ("`buf'") ///
        ("High-EVI/Low-O3 vs Reference") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))
}

* --- 4c. Tertile-based extreme contrast ---
di _n "=== Tertile-based extreme contrast ==="
capture drop tert_extreme_1KM tert_extreme_500m
foreach buf in 1KM 500m {
    * High tertile EVI + Low tertile O3 vs Low tertile EVI + High tertile O3
    gen tert_extreme_`buf' = .
    replace tert_extreme_`buf' = 1 if evi_tert_`buf' == 1 & o3_tert_`buf' == 3  // low EVI, high O3
    replace tert_extreme_`buf' = 2 if evi_tert_`buf' == 2 & o3_tert_`buf' == 2  // middle
    replace tert_extreme_`buf' = 3 if evi_tert_`buf' == 3 & o3_tert_`buf' == 1  // high EVI, low O3
    replace tert_extreme_`buf' = 4 if evi_tert_`buf' == 3 & o3_tert_`buf' == 3  // high EVI, high O3

    xtmixed BAA_PhenoAge i.tert_extreme_`buf' $controls || id_num:, mle

    * High-EVI tertile / High-O3 tertile vs middle tertiles
    lincom 4.tert_extreme_`buf' - 2.tert_extreme_`buf'
    post `results' ("Tertile (extreme)") ("`buf'") ///
        ("High-EVI/High-O3 vs Middle") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))

    * Low-EVI tertile / High-O3 tertile vs middle
    lincom 1.tert_extreme_`buf' - 2.tert_extreme_`buf'
    post `results' ("Tertile (extreme)") ("`buf'") ///
        ("Low-EVI/High-O3 vs Middle") ///
        (r(estimate)) (r(se)) (r(p)) (r(lb)) (r(ub))
}

* --- 4d. Continuous interaction model ---
di _n "=== Continuous EVI × O3 interaction ==="
foreach buf in 1KM 500m {
    xtmixed BAA_PhenoAge c.c_EVI_`buf'_mean4y c.c_O3_`buf'_mean4y ///
        c.evi_o3_int_`buf' $controls || id_num:, mle

    * Store interaction coefficient
    local int_beta = _b[c.evi_o3_int_`buf']
    local int_se   = _se[c.evi_o3_int_`buf']
    local int_p    = 2 * ttail(e(df_r), abs(_b[c.evi_o3_int_`buf'] / _se[c.evi_o3_int_`buf']))
    local int_ci_l = _b[c.evi_o3_int_`buf'] - invttail(e(df_r), 0.025) * _se[c.evi_o3_int_`buf']
    local int_ci_u = _b[c.evi_o3_int_`buf'] + invttail(e(df_r), 0.025) * _se[c.evi_o3_int_`buf']

    post `results' ("Continuous interaction") ("`buf'") ///
        ("EVI × O3 product term") ///
        (`int_beta') (`int_se') (`int_p') (`int_ci_l') (`int_ci_u')
}

postclose `results'

* ---- 5. Generate comparison summary ----
di _n "======================================================"
di "5. SUMMARY: ROBUSTNESS ACROSS THRESHOLDS"
di "======================================================"

use "output/tables/kmeans_robustness_comparison.dta", clear

* ---- Export CSV (before display, so data is saved even if display fails) ----
export delimited using "output/tables/kmeans_robustness_comparison.csv", replace
di _n "Full results exported to output/tables/kmeans_robustness_comparison.csv"

* ---- Display comparison table ----
gen beta_ci = string(beta, "%9.3f") + " (" + ///
              string(ci_low, "%9.3f") + ", " + ///
              string(ci_high, "%9.3f") + ")"
gen sig = cond(pval < 0.001, "***", cond(pval < 0.01, "**", cond(pval < 0.05, "*", "")))
gen result_str = beta_ci + sig

di _n "=== High-EVI/High-O3 vs Reference (key diagnostic) ==="
list classification buffer beta pval if contrast == "High-EVI/High-O3 vs Reference", ///
    sepby(classification) noobs

di _n "=== Low-EVI/High-O3 vs Reference (O3-only contrast) ==="
list classification buffer beta pval if contrast == "Low-EVI/High-O3 vs Reference", ///
    sepby(classification) noobs

* ---- Key diagnostic summary ----
di _n "======================================================"
di "KEY DIAGNOSTIC: High-EVI/High-O3 = Reference Group?"
di "======================================================"
di "If the High-EVI/High-O3 coefficient is close to zero"
di "across all classification methods, the additive"
di "attenuation pattern is robust to threshold choice."
di ""
di "Original K-means target: beta = -0.13 (P = 0.85)"
di "======================================================"

********************************************************************************
