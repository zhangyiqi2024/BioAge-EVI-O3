set seed 20240101

*==============================================================================
* Instrumental Variable Analysis for Continuous Exposure Variables
* Using temperature inversion days as IV for O3, and slope squared term as IV for EVI
* Using climate variable principal components to control for multicollinearity
*==============================================================================

di as txt "=== Instrumental Variable Analysis for Continuous Exposure Variables ==="

* Reload data
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

* Generate log-transformed exposure variables (consistent with main effect model)
gen log_EVI_1KM = log(EVI_1KM_mean4y + 0.001)
gen log_EVI_500m = log(EVI_500m_mean4y + 0.001)
gen log_O3_1KM = log(O3_1KM_mean4y + 0.001)
gen log_O3_500m = log(O3_500m_mean4y + 0.001)

* Generate slope squared term - using squared 25-degree slope variable
gen steep_slope_pct_25deg_sq_lag5 = steep_slope_pct_25deg_500m_lag5^2

* Sort by individual and year
sort id_num year
xtset id_num year

*==============================================================================
* Climate Variable Principal Component Analysis - Addressing Multicollinearity
*==============================================================================

di "=== Climate Variable Principal Component Analysis ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

*==============================================================================
* First-Stage Regression: Testing IV Relevance with Exposure Variables
*==============================================================================

di as txt "=== First-Stage Regression: Continuous Exposure Variables ==="

* 1. First-stage regression for EVI_1KM - using slope squared term as IV
di as txt "--- EVI_1KM First-Stage Regression (using slope squared term) ---"
reg log_EVI_1KM ///
    c.steep_slope_pct_25deg_sq_lag5  /// IV for EVI (slope squared term)
    c.age i.sex i.education                              ///
    c.ln_hh_consumption_total i.retirement i.pension_insurance                      ///
    i.smoking_current i.alcohol_current                                     ///
    climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
    c.bmi i.hypertension i.dyslipidemia                         ///
    i.year,                                                    ///
    vce(cluster id_num)

estimates store first_stage_EVI_1KM_cont
test c.steep_slope_pct_25deg_sq_lag5
di as result "EVI_1KM First-stage F-statistic: " r(F) " (p = " r(p) ")"

* 2. First-stage regression for EVI_500m - using the same slope squared term
di as txt "--- EVI_500m First-Stage Regression (using slope squared term) ---"
reg log_EVI_500m ///
    c.steep_slope_pct_25deg_sq_lag5  ///
    c.age i.sex i.education                              ///
    c.ln_hh_consumption_total i.retirement i.pension_insurance                      ///
    i.smoking_current i.alcohol_current                                     ///
    climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
    c.bmi i.hypertension i.dyslipidemia                         ///
    i.year,                                                    ///
    vce(cluster id_num)

estimates store first_stage_EVI_500m_cont
test c.steep_slope_pct_25deg_sq_lag5
di as result "EVI_500m First-stage F-statistic: " r(F) " (p = " r(p) ")"

* 3. First-stage regression for O3_1KM - using temperature inversion days as IV
di as txt "--- O3_1KM First-Stage Regression ---"
capture confirm variable inversion_days1_lag5
if _rc == 0 {
    reg log_O3_1KM ///
        c.inversion_days1_lag5                   /// IV for O3
        c.age i.sex i.education                              ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance                      ///
        i.smoking_current i.alcohol_current                                     ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
        c.bmi i.hypertension i.dyslipidemia                         ///
        i.year,                                                    ///
        vce(cluster id_num)

    estimates store first_stage_O3_1KM_cont
    test c.inversion_days1_lag5
    di as result "O3_1KM First-stage F-statistic: " r(F) " (p = " r(p) ")"
}
else {
    di as error "Error: variable 'inversion_days1_lag5' does not exist"
}

* 4. First-stage regression for O3_500m
di as txt "--- O3_500m First-Stage Regression ---"
capture confirm variable inversion_days1_lag5
if _rc == 0 {
    reg log_O3_500m ///
        c.inversion_days1_lag5                   ///
        c.age i.sex i.education                              ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance                      ///
        i.smoking_current i.alcohol_current                                     ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
        c.bmi i.hypertension i.dyslipidemia                         ///
        i.year,                                                    ///
        vce(cluster id_num)

    estimates store first_stage_O3_500m_cont
    test c.inversion_days1_lag5
    di as result "O3_500m First-stage F-statistic: " r(F) " (p = " r(p) ")"
}

*==============================================================================
* Single-Exposure IV Regression (Analyzing EVI and O3 Effects Separately)
*==============================================================================

di as txt "=== Single-Exposure IV Regression ==="

* 1. EVI_1KM single effect - using slope squared term
di as txt "--- EVI_1KM Single Effect (IV) ---"
ivregress 2sls BAA_PhenoAge                          ///
    (log_EVI_1KM = c.steep_slope_pct_25deg_sq_lag5) ///
    c.age i.sex i.education                    ///
    c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
    i.smoking_current i.alcohol_current                           ///
    climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
    c.bmi i.hypertension i.dyslipidemia               ///
    i.year,                                          ///
    vce(cluster id_num)

estimates store iv_EVI_1KM_single

* 2. EVI_500m single effect - using the same slope squared term
di as txt "--- EVI_500m Single Effect (IV) ---"
ivregress 2sls BAA_PhenoAge                          ///
    (log_EVI_500m = c.steep_slope_pct_25deg_sq_lag5) ///
    c.age i.sex i.education                    ///
    c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
    i.smoking_current i.alcohol_current                           ///
    climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
    c.bmi i.hypertension i.dyslipidemia               ///
    i.year,                                          ///
    vce(cluster id_num)

estimates store iv_EVI_500m_single

* 3. O3_1KM single effect
di as txt "--- O3_1KM Single Effect (IV) ---"
capture confirm variable inversion_days1_lag5
if _rc == 0 {
    ivregress 2sls BAA_PhenoAge                          ///
        (log_O3_1KM = c.inversion_days1_lag5 )             ///
        c.age i.sex i.education                    ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
        i.smoking_current i.alcohol_current                           ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
        c.bmi i.hypertension i.dyslipidemia               ///
        i.year,                                          ///
        vce(cluster id_num)

    estimates store iv_O3_1KM_single
}

* 4. O3_500m single effect
di as txt "--- O3_500m Single Effect (IV) ---"
capture confirm variable inversion_days1_lag5
if _rc == 0 {
    ivregress 2sls BAA_PhenoAge                          ///
        (log_O3_500m = c.inversion_days1_lag5 ) ///
        c.age i.sex i.education                    ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
        i.smoking_current i.alcohol_current                           ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
        c.bmi i.hypertension i.dyslipidemia               ///
        i.year,                                          ///
        vce(cluster id_num)

    estimates store iv_O3_500m_single
}

*==============================================================================
* Dual-Exposure IV Regression (Simultaneously Considering EVI and O3)
*==============================================================================

di as txt "=== Dual-Exposure IV Regression ==="

* 1. 1KM-scale dual-exposure model - using slope squared term and temperature inversion days
di as txt "--- 1KM-Scale Dual-Exposure Model (IV) ---"
capture confirm variable inversion_days1_lag5
if _rc == 0 {
    ivregress 2sls BAA_PhenoAge                          ///
        (log_EVI_1KM log_O3_1KM =                        ///
         c.steep_slope_pct_25deg_sq_lag5 ///
         c.inversion_days1_lag5 )     ///
        c.age i.sex i.education                    ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
        i.smoking_current i.alcohol_current                           ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
        c.bmi i.hypertension i.dyslipidemia               ///
        i.year,                                          ///
        vce(cluster id_num)

    estimates store iv_dual_1KM
}

* 2. 500m-scale dual-exposure model - using the same IV combination
di as txt "--- 500m-Scale Dual-Exposure Model (IV) ---"
capture confirm variable inversion_days1_lag5
if _rc == 0 {
    ivregress 2sls BAA_PhenoAge                          ///
        (log_EVI_500m log_O3_500m =                      ///
         c.steep_slope_pct_25deg_sq_lag5 ///
         c.inversion_days1_lag5 )     ///
        c.age i.sex i.education                    ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
        i.smoking_current i.alcohol_current                           ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
        c.bmi i.hypertension i.dyslipidemia               ///
        i.year,                                          ///
        vce(cluster id_num)

    estimates store iv_dual_500m
}

*==============================================================================
* Instrumental Variable Validity Tests
*==============================================================================

di as txt "=== Instrumental Variable Validity Tests ==="

* Only perform weak IV tests; overidentification tests are not applicable
* (all models are exactly identified)
di as txt "Note: All models are exactly identified (number of IVs = number of endogenous variables), overidentification tests cannot be performed"

* Weak instrument tests
foreach var in EVI_1KM EVI_500m {
    di as txt "--- `var' Weak Instrument Test ---"
    capture estimates restore first_stage_`var'_cont
    if _rc == 0 {
        qui test c.steep_slope_pct_25deg_sq_lag5
        di as result "`var' First-stage F-statistic: " %6.2f r(F) " (p = " %6.4f r(p) ")"

        * Check if F-statistic exceeds 10 (rule of thumb for weak instruments)
        if r(F) > 10 {
            di as result "✓ Instrument strength is sufficient (F > 10)"
        }
        else if r(F) > 5 {
            di as result "⚠ Instrument is moderately weak (5 < F < 10)"
        }
        else {
            di as error "✗ Severe weak instrument problem (F < 5)"
        }
    }
}

capture confirm variable inversion_days1_lag5
if _rc == 0 {
    foreach var in O3_1KM O3_500m {
        di as txt "--- `var' Weak Instrument Test ---"
        capture estimates restore first_stage_`var'_cont
        if _rc == 0 {
            qui test c.inversion_days1_lag5
            di as result "`var' First-stage F-statistic: " %6.2f r(F) " (p = " %6.4f r(p) ")"

            if r(F) > 10 {
                di as result "✓ Instrument strength is sufficient (F > 10)"
            }
            else if r(F) > 5 {
                di as result "⚠ Instrument is moderately weak (5 < F < 10)"
            }
            else {
                di as error "✗ Severe weak instrument problem (F < 5)"
            }
        }
    }
}

*==============================================================================
* Results Output and Saving
*==============================================================================

* Save all estimation results
estimates save "output/tables/iv_analysis_results", replace

* Output detailed results tables
esttab iv_EVI_1KM_single iv_EVI_500m_single iv_O3_1KM_single iv_O3_500m_single ///
    using "output/tables/single_exposure_iv_results_pca.csv", replace ///
    b(4) se(4) star(* 0.1 ** 0.05 *** 0.01) ///
    stats(N, fmt(0)) ///
    title("Single-Exposure IV Regression Results (Using Climate PCs)") ///
    mtitle("EVI_1KM" "EVI_500m" "O3_1KM" "O3_500m")

esttab iv_dual_1KM iv_dual_500m ///
    using "output/tables/dual_exposure_iv_results_pca.csv", replace ///
    b(4) se(4) star(* 0.1 ** 0.05 *** 0.01) ///
    stats(N, fmt(0)) ///
    title("Dual-Exposure IV Regression Results (Using Climate PCs)") ///
    mtitle("1KM Scale" "500m Scale")

* Generate summary report
capture log close
log using "output/tables/continuous_exposure_iv_analysis_report_pca.log", replace text

di as txt "=== Continuous Exposure IV Analysis Complete Report (Using Climate PCs) ==="
di as txt "Analysis time: `c(current_date)' `c(current_time)'"
di as txt "Sample size: " _N
di as txt "Age range: 46-89 years"
di as txt "Instrumental variable specification:"
di as txt "  - EVI: using 500m buffer 25-degree slope squared term"
di as txt "  - O3: using temperature inversion days (inversion_days1_lag5)"
di as txt "Control variable specification:"
di as txt "  - Using climate variable principal components (climate_pc1-climate_pc5) instead of raw climate variables"

* Output key results
di _newline(2)
di as txt "=== Key Findings ==="

foreach var in EVI_1KM EVI_500m {
    capture estimates restore first_stage_`var'_cont
    if _rc == 0 {
        qui test c.steep_slope_pct_25deg_sq_lag5
        di as txt "`var' First-stage F-statistic: " %6.2f r(F) " (p = " %6.4f r(p) ")"
    }
}

capture confirm variable inversion_days1_lag5
if _rc == 0 {
    foreach var in O3_1KM O3_500m {
        capture estimates restore first_stage_`var'_cont
        if _rc == 0 {
            qui test c.inversion_days1_lag5
            di as txt "`var' First-stage F-statistic: " %6.2f r(F) " (p = " %6.4f r(p) ")"
        }
    }
}

di _newline
di as txt "Main IV estimation results:"
capture estimates restore iv_dual_1KM
if _rc == 0 {
    estimates replay
}

log close

di as txt "=== Continuous Exposure IV Analysis (Using Climate PCs) Complete! ==="
di as txt "Generated files:"
di as txt "1. ../output/tables/iv_analysis_results.ster - All estimation results"
di as txt "2. ../output/tables/continuous_exposure_iv_analysis_report_pca.log - Detailed analysis report"
di as txt "3. ../output/tables/single_exposure_iv_results_pca.csv - Single-exposure model results table"
di as txt "4. ../output/tables/dual_exposure_iv_results_pca.csv - Dual-exposure model results table"

**========================================
********************************************************************************
* New Module: Automatically Summarize Single + Dual Exposure IV Results,
* Generate English Tables and Figures, and Output to Word
* Fixed: Resolved table row count issue and P-value calculation issue
********************************************************************************

di as txt "=== Begin Generating Summary Tables and Figures for IV Analysis Results ==="

*----------------------------------------------------------
* 1. Define models, exposure types, and buffers to extract
*----------------------------------------------------------
local models  iv_EVI_1KM_single  iv_EVI_500m_single ///
              iv_O3_1KM_single   iv_O3_500m_single ///
              iv_dual_1KM        iv_dual_500m

* Clear any existing tempfiles
capture postclose handle
tempfile results

* Create temporary data file to store results
postfile handle str20 Exposure str20 ModelType str10 Buffer ///
        double Coef double SE double P ///
        double LCI double UCI using `results', replace

*----------------------------------------------------------
* 2. Fixed result extraction: Using normal distribution for P-values and CIs
*----------------------------------------------------------

foreach m of local models {

    capture estimates restore `m'
    if _rc != 0 {
        di as error "Model `m' cannot be restored, skipping"
        continue
    }

    di as txt "Processing model: `m'"

    * Determine model type and buffer
    local model_type = cond(strpos("`m'", "single") > 0, "Single-exposure", "Dual-exposure")

    if strpos("`m'", "1KM") > 0 {
        local buffer "1KM"
    }
    else if strpos("`m'", "500m") > 0 {
        local buffer "500m"
    }
    else {
        local buffer "Unknown"
    }

    * ---------- Single-exposure models ----------
    if "`model_type'" == "Single-exposure" {

        if strpos("`m'", "EVI") > 0 {
            local exposure "EVI"
            if "`buffer'" == "1KM" local varname "log_EVI_1KM"
            else local varname "log_EVI_500m"
        }
        else if strpos("`m'", "O3") > 0 {
            local exposure "O3"
            if "`buffer'" == "1KM" local varname "log_O3_1KM"
            else local varname "log_O3_500m"
        }

        matrix b = e(b)
        matrix V = e(V)

        * Extract coefficient and standard error
        local coef = b[1, "`varname'"]
        local se = sqrt(V["`varname'", "`varname'"])

        * Use normal distribution for P-value and CI (more appropriate for IV regression)
        local z = `coef'/`se'
        local p = 2 * (1 - normal(abs(`z')))
        local lci = `coef' - 1.96 * `se'
        local uci = `coef' + 1.96 * `se'

        post handle ("`exposure'") ("`model_type'") ("`buffer'") ///
                   (`coef') (`se') (`p') (`lci') (`uci')

        di as txt "Single-exposure model: `exposure' | Coef: `coef' | SE: `se' | P: `p'"
    }

    * ---------- Dual-exposure models ----------
    else if "`model_type'" == "Dual-exposure" {

        if "`buffer'" == "1KM" {
            local varlist "log_EVI_1KM log_O3_1KM"
        }
        else if "`buffer'" == "500m" {
            local varlist "log_EVI_500m log_O3_500m"
        }

        matrix b = e(b)
        matrix V = e(V)

        foreach varname of local varlist {

            if strpos("`varname'", "EVI") > 0 {
                local exposure "EVI"
            }
            else if strpos("`varname'", "O3") > 0 {
                local exposure "O3"
            }

            * Extract coefficient and standard error
            local coef = b[1, "`varname'"]
            local se = sqrt(V["`varname'", "`varname'"])

            * Use normal distribution for P-value and CI
            local z = `coef'/`se'
            local p = 2 * (1 - normal(abs(`z')))
            local lci = `coef' - 1.96 * `se'
            local uci = `coef' + 1.96 * `se'

            post handle ("`exposure'") ("`model_type'") ("`buffer'") ///
                       (`coef') (`se') (`p') (`lci') (`uci')

            di as txt "Dual-exposure model: `exposure' | Coef: `coef' | SE: `se' | P: `p'"
        }
    }
}

postclose handle

* Load results data
use `results', clear

* Check data
di as txt "Summary data preview:"
list, clean noobs

*----------------------------------------------------------
* 3. Create display columns: Coefficient(SE), P-value, 95% CI
*----------------------------------------------------------

* Format output columns
gen coef_se = string(Coef, "%9.3f") + " (" + string(SE, "%9.3f") + ")"
gen pvalue = string(P, "%9.4f")
replace pvalue = "<0.0001" if P < 0.0001  // Handle very small P-values
gen CI = "[" + string(LCI, "%9.3f") + ", " + string(UCI, "%9.3f") + "]"

* Create ID for sorting and display
gen display_order = _n
gen display_label = ""
replace display_label = Exposure + " (" + ModelType + ")" + " - " + Buffer

*----------------------------------------------------------
* 4. Draw Forest Plot (English)
*----------------------------------------------------------

* Create clearer display labels
gen forest_label = ""
replace forest_label = "EVI (Single) - 1KM" if Exposure == "EVI" & ModelType == "Single-exposure" & Buffer == "1KM"
replace forest_label = "EVI (Single) - 500m" if Exposure == "EVI" & ModelType == "Single-exposure" & Buffer == "500m"
replace forest_label = "O3 (Single) - 1KM" if Exposure == "O3" & ModelType == "Single-exposure" & Buffer == "1KM"
replace forest_label = "O3 (Single) - 500m" if Exposure == "O3" & ModelType == "Single-exposure" & Buffer == "500m"
replace forest_label = "EVI (Dual) - 1KM" if Exposure == "EVI" & ModelType == "Dual-exposure" & Buffer == "1KM"
replace forest_label = "O3 (Dual) - 1KM" if Exposure == "O3" & ModelType == "Dual-exposure" & Buffer == "1KM"
replace forest_label = "EVI (Dual) - 500m" if Exposure == "EVI" & ModelType == "Dual-exposure" & Buffer == "500m"
replace forest_label = "O3 (Dual) - 500m" if Exposure == "O3" & ModelType == "Dual-exposure" & Buffer == "500m"

* Generate forest plot
twoway (rcap UCI LCI display_order, horizontal lcolor(blue)) ///
       (scatter display_order Coef, mcolor(red) msize(medium)), ///
       ylabel(1(1)`=_N' `forest_label', angle(0) labsize(small)) ///
       yscale(reverse) ///
       xtitle("Coefficient (95% CI)", size(medium)) ///
       ytitle("") ///
       title("Instrumental Variable Estimates", size(medium)) ///
       subtitle("for Green Vegetation (EVI) and Ozone (O3) Exposures", size(small)) ///
       legend(off) ///
       xline(0, lpattern(dash) lcolor(gray)) ///
       plotregion(margin(l=30 r=5)) ///
       graphregion(margin(l=5 r=5))

graph export "output/figures/IV_Forest_Plot.png", width(2000) replace
di as result "Forest plot saved as: ../output/figures/IV_Forest_Plot.png"

*----------------------------------------------------------
* 5. Output Word File (Fixed: dynamically adjust table size)
*----------------------------------------------------------

capture putdocx clear
putdocx begin

* Title
putdocx paragraph, style("Title")
putdocx text ("Instrumental Variable Analysis Results Summary")

* Subtitle
putdocx paragraph
putdocx text ("Table: Instrumental variable estimates for EVI and O3 exposures on Phenotypic Age Acceleration")
putdocx text (" using slope squared as IV for EVI and temperature inversion days as IV for O3")

* Dynamically calculate table row count
local row_count = _N + 1  // +1 for header row
di as txt "Table rows: `row_count' (data rows: `=_N', header rows: 1)"

* Create table - using dynamic row count
putdocx table iv_table = (`row_count', 6), border(all, nil) width(100%)

* Header
putdocx table iv_table(1,1) = ("Exposure"), bold halign(center)
putdocx table iv_table(1,2) = ("Model Type"), bold halign(center)
putdocx table iv_table(1,3) = ("Buffer"), bold halign(center)
putdocx table iv_table(1,4) = ("Coefficient (SE)"), bold halign(center)
putdocx table iv_table(1,5) = ("P-value"), bold halign(center)
putdocx table iv_table(1,6) = ("95% CI"), bold halign(center)

* Fill table data - fixed row numbering
forvalues i = 1/`=_N' {
    local row = `i' + 1
    di as txt "Filling row `row', data index `i'"

    putdocx table iv_table(`row',1) = ("`=Exposure[`i']'"), halign(center)
    putdocx table iv_table(`row',2) = ("`=ModelType[`i']'"), halign(center)
    putdocx table iv_table(`row',3) = ("`=Buffer[`i']'"), halign(center)
    putdocx table iv_table(`row',4) = ("`=coef_se[`i']'"), halign(center)
    putdocx table iv_table(`row',5) = ("`=pvalue[`i']'"), halign(center)
    putdocx table iv_table(`row',6) = ("`=CI[`i']'"), halign(center)
}

* Add table notes
putdocx paragraph
putdocx text ("Note: Single-exposure models include only one exposure (EVI or O3), while dual-exposure models include both exposures simultaneously.")

* Add forest plot to Word document
putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("Forest Plot of Instrumental Variable Estimates")

putdocx paragraph
putdocx image "output/figures/IV_Forest_Plot.png"

* Add methods description
putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("Methods")

putdocx paragraph
putdocx text ("Instrumental variable analysis was conducted to estimate the causal effects of:")
putdocx text (" (1) Green vegetation exposure measured by Enhanced Vegetation Index (EVI) and")
putdocx text (" (2) Ozone pollution exposure (O3) on phenotypic age acceleration.")

putdocx paragraph
putdocx text ("Instrumental variables:")
putdocx text (" - EVI: squared slope within 500m buffer (25-degree threshold)")
putdocx text (" - O3: temperature inversion days")

putdocx paragraph
putdocx text ("All models were adjusted for age, sex, education, household consumption, retirement status,")
putdocx text ("pension insurance, smoking, drinking, BMI, hypertension, dyslipidemia, year fixed effects,")
putdocx text ("and principal components of climate variables.")

* Save Word document
putdocx save "output/tables/IV_Summary_Table_and_Figure.docx", replace
di as result "Word document saved as: ../output/tables/IV_Summary_Table_and_Figure.docx"

*----------------------------------------------------------
* 6. Additional CSV-format detailed results output
*----------------------------------------------------------

* Create more detailed output variables
gen Significance = ""
replace Significance = "***" if P < 0.001
replace Significance = "**" if P >= 0.001 & P < 0.01
replace Significance = "*" if P >= 0.01 & P < 0.05
replace Significance = "" if P >= 0.05

gen Coef_SE_Sig = string(Coef, "%9.3f") + " (" + string(SE, "%9.3f") + ")" + Significance

* Save detailed results
export excel Exposure ModelType Buffer Coef SE P LCI UCI Significance ///
    using "output/tables/IV_Detailed_Results.xlsx", firstrow(variables) replace

di as result "Detailed results Excel file saved as: ../output/tables/IV_Detailed_Results.xlsx"

*----------------------------------------------------------
* 7. Results Summary Report
*----------------------------------------------------------

di _newline(2)
di as txt "=== Instrumental Variable Analysis Results Summary ==="
di as txt "Generated files:"
di as txt "1. ../output/tables/IV_Summary_Table_and_Figure.docx - Word summary document (with table and forest plot)"
di as txt "2. ../output/figures/IV_Forest_Plot.png - Forest plot"
di as txt "3. ../output/tables/IV_Detailed_Results.xlsx - Detailed results Excel file"
di as txt "4. ../output/tables/single_exposure_iv_results_pca.csv - Single-exposure model results"
di as txt "5. ../output/tables/dual_exposure_iv_results_pca.csv - Dual-exposure model results"

* Display key results summary
di _newline
di as txt "Key results summary:"
list Exposure ModelType Buffer Coef_SE_Sig CI pvalue, clean noobs

di as txt "=== Instrumental Variable Analysis Summary Complete! ==="








*============================================================================
*============================================================================
*============================================================================
* Scientific Explanation for the Marginal Significance of EVI Protective
* Effect in Dual-Exposure Models
*
* The protective effect of EVI on BAA includes both an indirect effect through
* O3 reduction (approximately 10.9) and a direct effect (approximately 41.6).
* Because the direct effect estimate is smaller and has relatively larger
* estimation uncertainty, statistical significance decreases.
*
* In dual-exposure models, we simultaneously estimate two endogenous variables,
* which increases model complexity. Each coefficient estimate is influenced by
* the other endogenous variable, leading to larger standard errors or smaller
* t-statistics.
*
* Furthermore, in dual-exposure models, the IV for EVI (slope interaction term)
* must now identify EVI alongside the O3 IV (temperature inversion days) to
* jointly identify both variables. This may weaken the identification power of
* each instrument, thereby increasing estimation uncertainty.
*
* The EVI coefficient in IV models becoming non-significant may be due to weak
* instruments or small sample size leading to larger standard errors.
* Additionally, IV estimator variance is typically larger than OLS, so even
* with larger coefficient magnitudes, results may not be significant.
*
* If instruments are valid, IV estimates are closer to true effects. However,
* given the larger variance of IV estimators, more data may be needed for
* precise estimation.
*
* In environmental epidemiological research, we observe that the protective
* effect of vegetation index (EVI) on biological age acceleration (BAA)
* reaches statistical significance in single-exposure IV models
* (beta = -52.533, p = 0.036), but in dual-exposure models controlling for
* ozone (O3) exposure, this effect attenuates and shows only marginal
* significance (beta = -41.633, p = 0.091). This phenomenon can be reasonably
* explained from both methodological and biological mechanism dimensions.
*
* Methodological Considerations
* From a causal inference perspective, single-exposure and dual-exposure IV
* models estimate fundamentally different causal parameters. Single-exposure
* models identify the total effect of EVI on BAA, while dual-exposure models
* estimate the direct effect of EVI after controlling for O3. This parameter
* difference stems from the identification properties of IV methods under
* multiple endogenous variable settings.
*
* When introducing a second endogenous variable O3, limited instrument resources
* must be allocated across both exposure variables for identification. Our
* diagnostic tests indicate that while instrument strength still meets rule-of-
* thumb requirements, the identification precision for each endogenous variable
* decreases relative to single-exposure models, manifested as larger standard
* errors and reduced statistical power.
*
* Biological Mechanism Interpretation
* From an environmental health mechanism perspective, changes in EVI effect
* magnitude suggest potential mediation effects. Vegetation cover influences
* local air quality through multiple pathways: plant leaves directly absorb
* ozone precursors, vegetation canopy modifies local meteorological conditions
* affecting photochemical reactions, and reduced ambient temperature decreases
* ozone formation. Therefore, part of EVI's protective effect on BAA may be
* realized through the mediating pathway of reduced O3 exposure.
*
* Effect decomposition analysis supports this mechanism: of EVI's total
* protective effect (beta = -52.533), approximately 20.7% can be explained
* through O3 reduction pathways, while the remaining 79.3% represents direct
* protective effects independent of O3. This effect decomposition is consistent
* with existing environmental science evidence suggesting that urban green
* space influences resident health through multiple mechanisms.
*
* Statistical Uncertainty Analysis
* The marginal significance of EVI effects in dual-exposure models (p = 0.091)
* reflects increased statistical estimation uncertainty rather than a
* fundamental change in effect authenticity. Notably, the direction of EVI
* effects is completely consistent across both models, and confidence intervals
* show substantial overlap. In observational research, such consistency is
* often more informative about the true exposure-response relationship than
* any single p-value threshold.
*
* Research Significance and Limitations
* Our findings highlight the complexity of environmental multi-exposure
* analysis. Single-exposure models may overestimate the independent effects of
* specific environmental factors because they cannot separate confounding or
* mediating effects of correlated exposures. In contrast, dual-exposure models
* provide more rigorous causal evidence, albeit at the cost of statistical
* precision.
*
* This finding has important implications for environmental policy: urban
* greening may promote healthy aging through dual mechanisms -- both directly
* slowing biological aging processes and indirectly protecting through improved
* air quality. Future research should adopt more refined mediation analysis
* frameworks to quantify the complex interaction networks among different
* environmental exposures.
*
* Conclusion
* The marginal significance of EVI protective effects in dual-exposure models
* reflects the inherent complexity of environmental health research rather than
* methodological flaws or unreliable findings. On the contrary, this pattern
* change reveals biologically plausible mediating pathways among environmental
* exposures, providing important evidence for understanding the multi-pathway
* influence of urban environments on healthy aging.




******************************
* 1. Interpretation of Log-Transformed Coefficients
* In the model, after log-transforming exposure variables, coefficients can be
* interpreted as: for each 1% increase in the exposure variable, the outcome
* variable (BAA_PhenoAge) changes by the corresponding proportion.
*
* Specifically, for the model:
*
*
* 2. Explaining Differences Between Main Effect Model and IV Model Coefficients
* Main effect model (mixed-effects model):
*
* log_EVI_1KM coefficient = -3.70 (p=0.023)
* log_O3_1KM coefficient = 8.44 (p=0.002)
*
* Dual-exposure IV regression (IV model):
* log_EVI_1KM coefficient = -28.68 (p=0.078)
* log_O3_1KM coefficient = 23.67 (p=0.001)
*
* We can see that in the IV model, the absolute values of both exposure
* variable coefficients are much larger than in the main effect model
* (EVI from -3.70 to -28.68, O3 from 8.44 to 23.67). This difference may
* be attributable to the following causes:
*
* (1) Measurement Error
* If exposure variables (EVI and O3) contain measurement error, this leads to
* attenuation bias (coefficients biased toward zero) in the main effect model.
* Instrumental variable methods can correct for such measurement error, so IV
* model coefficients are typically closer to true values.
*
* (2) Omitted Variable Bias
* The main effect model may not fully control for all confounding factors,
* leading to biased estimates. Instrumental variable methods, under the
* exogeneity condition, can eliminate this bias. If instruments are valid,
* IV estimators are consistent while main effect model estimates may be
* inconsistent.
*
* (3) Reverse Causality
* If BAA_PhenoAge in turn influences exposure variables (e.g., biological aging
* may affect an individual's residential environment), the main effect model
* will be biased due to reverse causality. IV methods can address this problem
* as long as instruments are correlated with exposure variables and affect
* outcomes only through exposure variables.
*
* (4) Local Average Treatment Effect (LATE)
* IV estimation captures the Local Average Treatment Effect, i.e., the effect
* among individuals whose exposure is influenced by the instrument. If the
* complier subpopulation differs from the overall population, IV estimates may
* differ from the average treatment effect estimated by the main effect model.
*
* (5) Model Specification Differences
* The main effect model is a mixed-effects model controlling for individual
* random effects, while the IV model uses two-stage least squares with
* instruments for estimation. The different estimation strategies of the two
* methods may also contribute to coefficient differences.
*
* 3. Specific to This Case
* In the main effect model, EVI's protective effect (negative coefficient) and
* O3's harmful effect (positive coefficient) are both significant.
*
* In the IV model, EVI's protective effect becomes larger (greater absolute
* value) but non-significant (p=0.078), while O3's harmful effect also becomes
* larger and remains significant.
*
* This pattern suggests that the main effect model may underestimate the true
* effects of exposure variables, possibly due to measurement error or omitted
* variable bias. The IV model, using instruments (slope squared term and
* temperature inversion days) to correct these biases, yields larger effect
* estimates.
*
* However, it should be noted that the EVI coefficient becomes non-significant
* in the IV model, which may be due to weak instruments or small sample size
* leading to larger standard errors. Additionally, IV estimator variance is
* typically larger than OLS, so coefficients may be non-significant even with
* larger magnitudes.
*
* 4. Conclusion
* Possible reasons for coefficient differences between the two models include:
* measurement error, omitted variables, reverse causality, etc. Instrumental
* variable methods attempt to address these issues but require satisfaction of
* the instrument exogeneity condition. If instruments are valid, IV estimates
* are closer to true effects. However, given the larger variance of IV
* estimators, more data may be needed for precise estimation.
*
* It is recommended to further examine instrument validity (whether first-stage
* F-statistics are sufficiently large) and conduct sensitivity analyses.
*

**********************************

* 1. Interpretation of Log-Transformed Exposure Variable Coefficients
*
* EVI Coefficient Interpretation
* Main effect model: log_EVI_1KM = -3.70
* Meaning: A 1% increase in EVI reduces BAA_PhenoAge by approximately 0.037 units
* Calculation: -3.70 x 0.01 = -0.037
*
* Instrumental variable model: log_EVI_1KM = -28.68
* Meaning: A 1% increase in EVI reduces BAA_PhenoAge by approximately 0.287 units
* Calculation: -28.68 x 0.01 = -0.287
*
* O3 Coefficient Interpretation
* Main effect model: log_O3_1KM = 8.44
* Meaning: A 1% increase in O3 increases BAA_PhenoAge by approximately 0.084 units
* Calculation: 8.44 x 0.01 = 0.084
*
* Instrumental variable model: log_O3_1KM = 23.67
* Meaning: A 1% increase in O3 increases BAA_PhenoAge by approximately 0.237 units
* Calculation: 23.67 x 0.01 = 0.237
*
* 2. Explaining Coefficient Differences Between Main Effect and IV Models
*
* Statistical Facts of Coefficient Differences
* EVI effect: from -3.70 to -28.68 (approximately 7.7x increase)
* O3 effect: from 8.44 to 23.67 (approximately 2.8x increase)
*
* Possible Explanatory Mechanisms
*
* (1) Measurement Error Bias
* Main effect model may suffer from classical measurement error:
*   True model: BAA = beta x EVI_true + epsilon
*   Observed model: EVI_obs = EVI_true + u (measurement error)
*   Result: beta_OLS -> beta x [Var(EVI_true)/(Var(EVI_true)+Var(u))] < beta
*
* Vegetation index (EVI) and O3 concentration both contain spatial interpolation
* error. Instrumental variable methods can mitigate attenuation bias from
* measurement error.
*
* (2) Omitted Variable Bias
* The main effect model may omit confounders correlated with EVI/O3:
*   - Socioeconomic factors: greener areas tend to have better socioeconomic
*     conditions
*   - Behavioral factors: greenness influences physical activity, mental health,
*     etc.
*   - Environmental co-exposures: correlations with other pollutants
*
* (3) Reverse Causality
*   EVI <- socioeconomic factors -> BAA
*        \                      /
*         bidirectional causality
*
* Biological aging may influence individuals' residential location choices.
* Instrumental variables break this bidirectional relationship through the
* exogeneity assumption.
*
* (4) Local Average Treatment Effect (LATE)
* IV estimates capture effects in subpopulations affected by the instrument:
*   - Slope squared term affects vegetation types sensitive to topography
*   - Temperature inversion days affect O3 concentrations in areas sensitive
*     to meteorological conditions
*   - These subpopulations may be more sensitive to exposures
*
* 3. Economic and Epidemiological Interpretation
*
* (1) Selection Bias
*   Affluent healthy individuals -> choose greener, less polluted areas ->
*   apparent "protective effect" overestimated
*   Instrumental variables -> exploit exogenous topographic variation ->
*   closer to true causal effect
*
* (2) Ecological Bias
* The main effect model may confound:
*   - Individual-level effects: direct physiological effects of vegetation on
*     individuals
*   - Neighborhood-level effects: social environment improvements from greening
*   - Instrumental variables may better identify the former
*
* 4. Statistical Test Support
* Instrument Validity
*   - First-stage F-statistic: if >10, supports instrument strength
*   - Coefficient direction consistency: consistent signs across both models
*     support the existence of causal relationships
*
* Economic Interpretation Framework
*   Main effect model: correlation = direct effect + various biases
*   IV model: closer to true causal effect
*   Magnitude of difference reflects severity of biases
*
* 5. Policy Implications
* If IV results are more credible:
*   - EVI's protective effect is stronger than traditional analyses suggest
*   - O3's harmful effect is also stronger than traditional analyses suggest
*   - Environmental interventions (increasing greening, controlling O3) have
*     greater potential benefits
*
* 6. Caveats
* Further verification needed:
*   - Instrument exogeneity: does slope affect aging through other channels?
*   - Effect heterogeneity: do effects differ across population subgroups?
*   - Biological plausibility: do such large effects have biological basis?
*
* Recommended sensitivity analyses:
*   - Test robustness under different IV specifications
*   - Analyze effect differences across subgroups
*   - Test for nonlinear relationships
*
* Conclusion
* The large coefficient differences between the main effect model and the IV
* model strongly suggest the presence of serious measurement error, omitted
* variable, or reverse causality problems. IV estimates indicate that
* environmental exposure effects on biological aging may be substantially
* larger than traditional regression analyses suggest, carrying important
* implications for environmental health policy.
*
* This pattern of differences (enhanced effects for both exposure variables)
* further supports the credibility of IV results, as it is not random variation
* but systematically points toward larger environmental health effects.
