* =============================================================================
* Step 1: Main Analysis - Six Models (Panel Cross-Section: 2011→2015 delta_BAA_PhenoAge)
* =============================================================================

set seed 20240101

di as txt "=== Starting Main Analysis - Six Models (Panel Cross-Sectional Analysis) ==="

use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

sort id_num year

* PCA for climate variables
di "=== Climate Variable Principal Component Analysis ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score


* =============================================================================
* Model 1: EVI
* =============================================================================

* Model 1-1: 1KM EVI
di as txt "=== Model 1-1: EVI_1KM delta_BAA_PhenoAge ==="
preserve
keep if year == 2011

regress delta_BAA_PhenoAge i.traj15_EVI_1KM ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current i.hobbies ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.BA_PhenoAge c.bmi i.hypertension i.dyslipidemia ///
     [pweight=ipw_EVI_1KM], vce(robust)

estimates store model1_EVI_1KM
restore


* Model 1-2: 500m EVI
di as txt "=== Model 1-2: EVI_500m delta_BAA_PhenoAge ==="
preserve
keep if year == 2011

regress delta_BAA_PhenoAge i.traj15_EVI_500m ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current i.hobbies ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.BA_PhenoAge c.bmi i.hypertension i.dyslipidemia ///
     [pweight=ipw_EVI_500m], vce(robust)

estimates store model1_EVI_500m
restore


* =============================================================================
* Model 2: O3
* =============================================================================

* Model 2-1: 1KM O3
di as txt "=== Model 2-1: O3_1KM delta_BAA_PhenoAge ==="
preserve
keep if year == 2011

regress delta_BAA_PhenoAge i.traj15_O3_1KM ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current i.hobbies ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.BA_PhenoAge c.bmi i.hypertension i.dyslipidemia ///
     [pweight=ipw_O3_1KM], vce(robust)

estimates store model2_O3_1KM
restore


* Model 2-2: 500m O3
di as txt "=== Model 2-2: O3_500m delta_BAA_PhenoAge ==="
preserve
keep if year == 2011

regress delta_BAA_PhenoAge i.traj15_O3_500m ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current i.hobbies ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.BA_PhenoAge c.bmi i.hypertension i.dyslipidemia ///
     [pweight=ipw_O3_500m], vce(robust)

estimates store model2_O3_500m
restore


* =============================================================================
* Model 3: Joint Model EVI + O3
* =============================================================================

* Model 3-1: 1KM Joint
di as txt "=== Model 3-1: 1KM Joint Model ==="
preserve
keep if year == 2011

regress delta_BAA_PhenoAge i.traj15_EVI_1KM i.traj15_O3_1KM ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current i.hobbies ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.BA_PhenoAge c.bmi i.hypertension i.dyslipidemia ///
     [pweight=ipw_joint_1KM], vce(robust)

estimates store model3_joint_1KM
restore


* Model 3-2: 500m Joint
di as txt "=== Model 3-2: 500m Joint Model ==="
preserve
keep if year == 2011

regress delta_BAA_PhenoAge i.traj15_EVI_500m i.traj15_O3_500m ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current i.hobbies ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.BA_PhenoAge c.bmi i.hypertension i.dyslipidemia ///
     [pweight=ipw_joint_500m], vce(robust)

estimates store model3_joint_500m
restore



* =============================================================================
* Step 2: Extract Eight Exposure Coefficients, Format as English-Language Research Table
* Exposure labels changed to:
*   EVI(Single-exposure), O3(Single-exposure),
*   EVI(Dual-exposure),   O3(Dual-exposure)
* (Model column removed)
* =============================================================================

capture postutil clear
tempname mem

postfile `mem' ///
    str30 Exposure str10 Buffer ///
    double Coef double SE double P ///
    str30 CI ///
    using results_summary.dta, replace


* --------------------------------------------------------------------------
* 1) model1_EVI_1KM, exposure = 1.traj15_EVI_1KM → EVI(Single-exposure), 1 km
* --------------------------------------------------------------------------
estimates restore model1_EVI_1KM
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_EVI_1KM"]
scalar se   = sqrt(V["1.traj15_EVI_1KM","1.traj15_EVI_1KM"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str1 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("EVI(Single-exposure)") ("1 km") (coef) (se) (p) ("`ci_str1'")


* --------------------------------------------------------------------------
* 2) model1_EVI_500m, exposure = 1.traj15_EVI_500m → EVI(Single-exposure), 500 m
* --------------------------------------------------------------------------
estimates restore model1_EVI_500m
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_EVI_500m"]
scalar se   = sqrt(V["1.traj15_EVI_500m","1.traj15_EVI_500m"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str2 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("EVI(Single-exposure)") ("500 m") (coef) (se) (p) ("`ci_str2'")


* --------------------------------------------------------------------------
* 3) model2_O3_1KM, exposure = 1.traj15_O3_1KM → O3(Single-exposure), 1 km
* --------------------------------------------------------------------------
estimates restore model2_O3_1KM
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_O3_1KM"]
scalar se   = sqrt(V["1.traj15_O3_1KM","1.traj15_O3_1KM"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str3 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("O3(Single-exposure)") ("1 km") (coef) (se) (p) ("`ci_str3'")


* --------------------------------------------------------------------------
* 4) model2_O3_500m, exposure = 1.traj15_O3_500m → O3(Single-exposure), 500 m
* --------------------------------------------------------------------------
estimates restore model2_O3_500m
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_O3_500m"]
scalar se   = sqrt(V["1.traj15_O3_500m","1.traj15_O3_500m"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str4 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("O3(Single-exposure)") ("500 m") (coef) (se) (p) ("`ci_str4'")


* --------------------------------------------------------------------------
* 5) model3_joint_1KM, exposure = 1.traj15_EVI_1KM → EVI(Dual-exposure), 1 km
* --------------------------------------------------------------------------
estimates restore model3_joint_1KM
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_EVI_1KM"]
scalar se   = sqrt(V["1.traj15_EVI_1KM","1.traj15_EVI_1KM"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str5 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("EVI(Dual-exposure)") ("1 km") (coef) (se) (p) ("`ci_str5'")


* --------------------------------------------------------------------------
* 6) model3_joint_1KM, exposure = 1.traj15_O3_1KM → O3(Dual-exposure), 1 km
* --------------------------------------------------------------------------
estimates restore model3_joint_1KM
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_O3_1KM"]
scalar se   = sqrt(V["1.traj15_O3_1KM","1.traj15_O3_1KM"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str6 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("O3(Dual-exposure)") ("1 km") (coef) (se) (p) ("`ci_str6'")


* --------------------------------------------------------------------------
* 7) model3_joint_500m, exposure = 1.traj15_EVI_500m → EVI(Dual-exposure), 500 m
* --------------------------------------------------------------------------
estimates restore model3_joint_500m
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_EVI_500m"]
scalar se   = sqrt(V["1.traj15_EVI_500m","1.traj15_EVI_500m"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str7 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("EVI(Dual-exposure)") ("500 m") (coef) (se) (p) ("`ci_str7'")


* --------------------------------------------------------------------------
* 8) model3_joint_500m, exposure = 1.traj15_O3_500m → O3(Dual-exposure), 500 m
* --------------------------------------------------------------------------
estimates restore model3_joint_500m
matrix b = e(b)
matrix V = e(V)

scalar coef = b[1,"1.traj15_O3_500m"]
scalar se   = sqrt(V["1.traj15_O3_500m","1.traj15_O3_500m"])
scalar df   = e(df_r)
scalar tval = coef/se
scalar p    = 2*ttail(df, abs(tval))
scalar lci  = coef - 1.96*se
scalar uci  = coef + 1.96*se
local ci_str8 = "[" + string(lci,"%6.3f") + ", " + string(uci,"%6.3f") + "]"

post `mem' ("O3(Dual-exposure)") ("500 m") (coef) (se) (p) ("`ci_str8'")


postclose `mem'


* -----------------------------------------------------------------------------
* Step 3: Format Numeric Values to 3 Decimal Places, Export as Word Table
* -----------------------------------------------------------------------------
use results_summary.dta, clear
order Exposure Buffer Coef SE P CI

* Set numeric variables to display 3 decimal places
format Coef SE P %9.3f

* (Optional) Preview the table
list, sep(0)

* Use putdocx to export Word table
putdocx clear
putdocx begin

* Title
putdocx paragraph, style(Title)
putdocx text ("Associations of Green Space and Ozone Exposure with Phenotypic Aging Acceleration (Delta BAA)")

* Brief description
putdocx paragraph
putdocx text ("Table shows regression coefficients, standard errors, P-values, and 95% confidence intervals for single- and dual-exposure models at 1 km and 500 m buffers.")

* Generate table (include variable names as headers)
putdocx table tbl1 = data(Exposure Buffer Coef SE P CI), varnames

* Save Word file
putdocx save "output/tables/results_summary_DeltaBAA.docx", replace

di as txt "=== English Word table exported: results_summary_DeltaBAA.docx ==="
