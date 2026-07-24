* ----------------------------------------------------------------------------- 
* Nonlinearity Tests - Four Main Models (Continuous Exposure + Climate PCs)
* -----------------------------------------------------------------------------

set seed 20240101

di as txt "=== Begin Nonlinearity Tests - Four Main Models (Continuous Exposure, Climate PCs) ==="

* Reload data (working directory must be the project root: the folder
* containing CODE/, DATA/ and output/ - see CODE/00_master_run.do)
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

* -----------------------------------------------------------------------------
* Climate PCA - Address Multicollinearity
* -----------------------------------------------------------------------------
di "=== Climate PCA ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* Sort by individual and year
sort id_num year

* Standardize continuous variables (reduce multicollinearity)
foreach var in EVI_1KM_lag1 EVI_500m_lag1 O3_1KM_lag1 O3_500m_lag1 {
    summarize `var'
    gen `var'_std = (`var' - r(mean)) / r(sd)
    gen `var'_std_sq = `var'_std^2
    gen `var'_std_cu = `var'_std^3
}

***********************Model 1: EVI_1KM BAA_PhenoAge (Mixed-effects model)
di as txt "=== Model 1: EVI_1KM BAA_PhenoAge (Continuous + nonlinear terms, using climate PCs) ==="

mixed BAA_PhenoAge c.EVI_1KM_lag1_std c.EVI_1KM_lag1_std_sq c.EVI_1KM_lag1_std_cu ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///    
     i.year ///
     || id_num:, ///
     pweight(ipw_EVI_1KM) vce(robust)

estimates store model1_nonlinear_pca

* Nonlinearity test
test c.EVI_1KM_lag1_std_sq c.EVI_1KM_lag1_std_cu
local nonlinear_p_1 = r(p)
di as txt "EVI_1KMNonlinearity testP值: " as res %8.4f `nonlinear_p_1'

* Calculate AIC and BIC
estat ic
matrix ic_1 = r(S)
local aic_1 = ic_1[1,5]
local bic_1 = ic_1[1,6]

***********************Model 2: EVI_500m BAA_PhenoAge (Mixed-effects model)
di as txt "=== Model 2: EVI_500m BAA_PhenoAge (Continuous + nonlinear terms, using climate PCs) ==="

mixed BAA_PhenoAge c.EVI_500m_lag1_std c.EVI_500m_lag1_std_sq c.EVI_500m_lag1_std_cu ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///    
     i.year ///
     || id_num:, ///
     pweight(ipw_EVI_500m) vce(robust)

estimates store model2_nonlinear_pca

* Nonlinearity test
test c.EVI_500m_lag1_std_sq c.EVI_500m_lag1_std_cu
local nonlinear_p_2 = r(p)
di as txt "EVI_500mNonlinearity testP值: " as res %8.4f `nonlinear_p_2'

* Calculate AIC and BIC
estat ic
matrix ic_2 = r(S)
local aic_2 = ic_2[1,5]
local bic_2 = ic_2[1,6]

***********************Model 3: O3_1KM BAA_PhenoAge (Mixed-effects model)
di as txt "=== Model 3: O3_1KM BAA_PhenoAge (Continuous + nonlinear terms, using climate PCs) ==="

mixed BAA_PhenoAge c.O3_1KM_lag1_std c.O3_1KM_lag1_std_sq c.O3_1KM_lag1_std_cu ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///    
     i.year ///
     || id_num:, ///
     pweight(ipw_O3_1KM) vce(robust)

estimates store model3_nonlinear_pca

* Nonlinearity test
test c.O3_1KM_lag1_std_sq c.O3_1KM_lag1_std_cu
local nonlinear_p_3 = r(p)
di as txt "O3_1KMNonlinearity testP值: " as res %8.4f `nonlinear_p_3'

* Calculate AIC and BIC
estat ic
matrix ic_3 = r(S)
local aic_3 = ic_3[1,5]
local bic_3 = ic_3[1,6]

***********************Model 4: O3_500m BAA_PhenoAge (Mixed-effects model)
di as txt "=== Model 4: O3_500m BAA_PhenoAge (Continuous + nonlinear terms, using climate PCs) ==="

mixed BAA_PhenoAge c.O3_500m_lag1_std c.O3_500m_lag1_std_sq c.O3_500m_lag1_std_cu ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///    
     i.year ///
     || id_num:, ///
     pweight(ipw_O3_500m) vce(robust)

estimates store model4_nonlinear_pca

* Nonlinearity test
test c.O3_500m_lag1_std_sq c.O3_500m_lag1_std_cu
local nonlinear_p_4 = r(p)
di as txt "O3_500mNonlinearity testP值: " as res %8.4f `nonlinear_p_4'

* Calculate AIC and BIC
estat ic
matrix ic_4 = r(S)
local aic_4 = ic_4[1,5]
local bic_4 = ic_4[1,6]

* -----------------------------------------------------------------------------
* Nonlinearity test结果汇总 - Professional scientific paper format
* -----------------------------------------------------------------------------

di as txt "=== Nonlinearity test结果汇总（Using climate PCs）==="

putdocx clear
putdocx begin

putdocx paragraph, style("Heading1")
putdocx text ("Nonlinearity Assessment of Environmental Exposures on Biological Aging")

putdocx paragraph
putdocx text ("Table 1. Nonlinearity Test Results for Environmental Exposure Variables and Biological Age Acceleration"), bold

putdocx table tbl1 = (5,5), border(all, single) width(100%)
putdocx table tbl1(1,1) = ("Exposure"), bold halign(center)
putdocx table tbl1(1,2) = ("Buffer"), bold halign(center)
putdocx table tbl1(1,3) = ("Nonlinearity P-value"), bold halign(center)
putdocx table tbl1(1,4) = ("AIC"), bold halign(center)
putdocx table tbl1(1,5) = ("BIC"), bold halign(center)

putdocx table tbl1(2,1) = ("EVI"), halign(center)
putdocx table tbl1(2,2) = ("1 km"), halign(center)
putdocx table tbl1(2,3) = (strofreal(`nonlinear_p_1', "%8.4f")), halign(center)
putdocx table tbl1(2,4) = (strofreal(`aic_1', "%8.1f")), halign(center)
putdocx table tbl1(2,5) = (strofreal(`bic_1', "%8.1f")), halign(center)

putdocx table tbl1(3,1) = ("EVI"), halign(center)
putdocx table tbl1(3,2) = ("500 m"), halign(center)
putdocx table tbl1(3,3) = (strofreal(`nonlinear_p_2', "%8.4f")), halign(center)
putdocx table tbl1(3,4) = (strofreal(`aic_2', "%8.1f")), halign(center)
putdocx table tbl1(3,5) = (strofreal(`bic_2', "%8.1f")), halign(center)

putdocx table tbl1(4,1) = ("Ozone"), halign(center)
putdocx table tbl1(4,2) = ("1 km"), halign(center)
putdocx table tbl1(4,3) = (strofreal(`nonlinear_p_3', "%8.4f")), halign(center)
putdocx table tbl1(4,4) = (strofreal(`aic_3', "%8.1f")), halign(center)
putdocx table tbl1(4,5) = (strofreal(`bic_3', "%8.1f")), halign(center)

putdocx table tbl1(5,1) = ("Ozone"), halign(center)
putdocx table tbl1(5,2) = ("500 m"), halign(center)
putdocx table tbl1(5,3) = (strofreal(`nonlinear_p_4', "%8.4f")), halign(center)
putdocx table tbl1(5,4) = (strofreal(`aic_4', "%8.1f")), halign(center)
putdocx table tbl1(5,5) = (strofreal(`bic_4', "%8.1f")), halign(center)

* Add statistical notes
putdocx paragraph
putdocx text ("Notes: All models adjusted for age, sex, education, socioeconomic status, lifestyle factors, health conditions, and climate principal components. Mixed-effects models with individual random intercepts were used. P-values > 0.05 indicate no significant nonlinear relationship."), italic

putdocx save "output/tables/Nonlinearity_Test_Results_Formatted.docx", replace

* -----------------------------------------------------------------------------
* Generate nonlinear relationship visualization - Professional scientific paper format
* -----------------------------------------------------------------------------

di as txt "=== Generate nonlinear relationship visualization（Professional scientific paper format）==="

* EVINonlinear relationship plot
mixed BAA_PhenoAge c.EVI_1KM_lag1_std c.EVI_1KM_lag1_std_sq c.EVI_1KM_lag1_std_cu ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///    
     i.year ///
     || id_num:, ///
     pweight(ipw_EVI_1KM) vce(robust)

margins, at(EVI_1KM_lag1_std=(-2(0.5)2)) saving(margins_evi_pca, replace)

marginsplot, ///
    title("Nonlinear Association between Enhanced Vegetation Index and Biological Age Acceleration", size(medium)) ///
    ytitle("Predicted Phenotypic Age Acceleration (Years)", size(small)) ///
    xtitle("Standardized EVI (1 km buffer)", size(small)) ///
    recast(line) recastci(rarea) ///
    plot1opts(lcolor("0 100 0") lwidth(0.8)) ///
    ci1opts(fcolor("0 100 0%30")) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    xlabel(-2(1)2, labsize(small)) ///
    ylabel(, labsize(small)) ///
    legend(off)

graph export "output/figures/Fig1_EVI_Nonlinear_Relationship.tif", replace width(2000)

* O3Nonlinear relationship plot
mixed BAA_PhenoAge c.O3_1KM_lag1_std c.O3_1KM_lag1_std_sq c.O3_1KM_lag1_std_cu ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///    
     i.year ///
     || id_num:, ///
     pweight(ipw_O3_1KM) vce(robust)

margins, at(O3_1KM_lag1_std=(-2(0.5)2)) saving(margins_o3_pca, replace)

marginsplot, ///
    title("Nonlinear Association between Ozone Exposure and Biological Age Acceleration", size(medium)) ///
    ytitle("Predicted Phenotypic Age Acceleration (Years)", size(small)) ///
    xtitle("Standardized O{subscript:3} (1 km buffer)", size(small)) ///
    recast(line) recastci(rarea) ///
    plot1opts(lcolor("200 0 0") lwidth(0.8)) ///
    ci1opts(fcolor("200 0 0%30")) ///
    graphregion(color(white)) ///
    plotregion(color(white)) ///
    xlabel(-2(1)2, labsize(small)) ///
    ylabel(, labsize(small)) ///
    legend(off)

graph export "output/figures/Fig2_O3_Nonlinear_Relationship.tif", replace width(2000)

di as txt "=== Nonlinearity test完成 ==="
di as txt "Generated files:"
di as txt "1. Nonlinearity_Test_Results_Formatted.docx - Formatted test results table"
di as txt "2. Fig1_EVI_Nonlinear_Relationship.tif - EVINonlinear relationship plot"
di as txt "3. Fig2_O3_Nonlinear_Relationship.tif - O3Nonlinear relationship plot"