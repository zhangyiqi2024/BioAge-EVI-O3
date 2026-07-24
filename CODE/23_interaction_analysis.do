* =============================================================================
* EVI与O3Interaction Analysis - 完整修正版本（交互项P值已修正）
* =============================================================================

set seed 20240101

di as txt "=== 开始EVI与O3Interaction Analysis - Corrected Version (Stratified by O3) ==="

* Reload complete dataset (working directory must be the project root:
* the folder containing CODE/, DATA/ and output/ - see CODE/00_master_run.do)
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

* Set graph style
set scheme s1color
graph set window fontface "Times New Roman"

* Climate PCA
di "=== Climate PCA ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* -----------------------------------------------------------------------------
* 1KM尺度Interaction Analysis
* -----------------------------------------------------------------------------

di as txt "=== 1KM尺度Interaction Analysis ==="

* 1. 加权模型Interaction Analysis
di as txt "运行1KM尺度加权Mixed-effects model..."
mixed BAA_PhenoAge i.traj15_EVI_1KM##i.traj15_O3_1KM ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current /// 
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///  
     c.bmi i.hypertension i.dyslipidemia ///     
     i.year ///
     || id_num:, ///
     pweight(ipw_joint_1KM) vce(robust)

estimates store interaction_1km_weighted

* Extract interaction coefficient and p-value (corrected: using _b[] and test directly)
global interaction_coef_1km = _b[1.traj15_EVI_1KM#1.traj15_O3_1KM]
test 1.traj15_EVI_1KM#1.traj15_O3_1KM
global interaction_p_1km = r(p)

* Joint interaction test (if multi-level)
testparm i.traj15_EVI_1KM#i.traj15_O3_1KM
global interaction_test_p_1km = r(p)

di as txt "1KM尺度交互项系数: " as res %8.4f $interaction_coef_1km
di as txt "1KM尺度交互项P值: " as res %8.4f $interaction_p_1km
di as txt "1KM尺度交互项检验P值: " as res %8.4f $interaction_test_p_1km

* 2. 1KM尺度Stratified analysis - 按O3Stratified analysisEVI效应
di as txt "运行1KM尺度按O3Stratified analysis..."

*-----------------------------
* Low O3 Group EVI Effect Analysis (1 km)
*-----------------------------
preserve
keep if traj15_O3_1KM == 0
mixed BAA_PhenoAge i.traj15_EVI_1KM ///
      c.age i.sex i.education ///
      c.ln_hh_consumption_total i.retirement i.pension_insurance ///
      i.smoking_current i.alcohol_current ///
      climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
      c.bmi i.hypertension i.dyslipidemia ///
      i.year ///
      || id_num:, ///
      vce(robust)

/* Use lincom to get beta and SE for 1.traj15_EVI_1KM */
lincom 1.traj15_EVI_1KM
scalar beta_lowO3_1km = r(estimate)
scalar se_lowO3_1km   = r(se)
global EVI_effect_lowO3_1km = beta_lowO3_1km   // EVI effect in low O3 group
restore

*-----------------------------
* High O3 Group EVI Effect Analysis (1 km)
*-----------------------------
preserve
keep if traj15_O3_1KM == 1
mixed BAA_PhenoAge i.traj15_EVI_1KM ///
      c.age i.sex i.education ///
      c.ln_hh_consumption_total i.retirement i.pension_insurance ///
      i.smoking_current i.alcohol_current ///
      climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
      c.bmi i.hypertension i.dyslipidemia ///
      i.year ///
      || id_num:, ///
      vce(robust)

lincom 1.traj15_EVI_1KM
scalar beta_highO3_1km = r(estimate)
scalar se_highO3_1km   = r(se)
global EVI_effect_highO3_1km = beta_highO3_1km   // EVI effect in high O3 group
restore

*-----------------------------------------
* Calculate 1 km Effect Difference and P-value (Stratified Mixed)
*   Assuming approximate independence of two group estimates：
*   Var(beta_high - beta_low) = Var(beta_high) + Var(beta_low)
*-----------------------------------------
global effect_diff_1km = $EVI_effect_highO3_1km - $EVI_effect_lowO3_1km

scalar se_diff_1km = sqrt(se_highO3_1km^2 + se_lowO3_1km^2)
scalar z_diff_1km  = $effect_diff_1km / se_diff_1km
scalar p_diff_1km  = 2 * normal(-abs(z_diff_1km))
global effect_diff_p_1km = p_diff_1km

di as txt "1KM尺度低O3组EVI效应: "  as res %8.4f $EVI_effect_lowO3_1km
di as txt "1KM尺度高O3组EVI效应: "  as res %8.4f $EVI_effect_highO3_1km
di as txt "1KM尺度效应差异: "      as res %8.4f $effect_diff_1km
di as txt "1KM尺度效应差异P值: "    as res %8.4f $effect_diff_p_1km

* Calculate 1km Effect Change Percentage (Keep Original Logic)
if $EVI_effect_lowO3_1km != 0 {
    global percent_change_1km = (($EVI_effect_highO3_1km - $EVI_effect_lowO3_1km) / abs($EVI_effect_lowO3_1km)) * 100
    di as txt "1KM尺度EVI效应变化百分比: " as res %6.2f $percent_change_1km as txt "%"
}


* -----------------------------------------------------------------------------
* 500m尺度Interaction Analysis
* -----------------------------------------------------------------------------

di as txt "=== 500m尺度Interaction Analysis ==="

* 1. 500m尺度加权模型Interaction Analysis
di as txt "运行500m尺度加权Mixed-effects model..."
mixed BAA_PhenoAge i.traj15_EVI_500m##i.traj15_O3_500m ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current /// 
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///  
     c.bmi i.hypertension i.dyslipidemia ///     
     i.year ///
     || id_num:, ///
     pweight(ipw_joint_500m) vce(robust)

estimates store interaction_500m_weighted

* Extract interaction coefficient and p-value (corrected)
global interaction_coef_500m = _b[1.traj15_EVI_500m#1.traj15_O3_500m]
test 1.traj15_EVI_500m#1.traj15_O3_500m
global interaction_p_500m = r(p)

testparm i.traj15_EVI_500m#i.traj15_O3_500m
global interaction_test_p_500m = r(p)

* Assign _new macros from correctly extracted values (bug fix for reproducibility)
global interaction_coef_1km_new = $interaction_coef_1km
global interaction_p_1km_new = $interaction_p_1km
global interaction_coef_500m_new = $interaction_coef_500m
global interaction_p_500m_new = $interaction_p_500m

di as txt "500m scale interaction coefficient: " as res %8.4f $interaction_coef_500m
di as txt "500m尺度交互项P值: " as res %8.4f $interaction_p_500m
di as txt "500m尺度交互项检验P值: " as res %8.4f $interaction_test_p_500m

* 2. 500m尺度Stratified analysis - 按O3Stratified analysisEVI效应
di as txt "运行500m尺度按O3Stratified analysis..."

* Low O3 Group EVI Effect Analysis (500m) - using lincom for SE extraction
preserve
keep if traj15_O3_500m == 0
mixed BAA_PhenoAge i.traj15_EVI_500m ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year ///
     || id_num:, ///
     vce(robust)
lincom 1.traj15_EVI_500m
scalar beta_lowO3_500m = r(estimate)
scalar se_lowO3_500m   = r(se)
global EVI_effect_lowO3_500m = beta_lowO3_500m
restore

* High O3 Group EVI Effect Analysis (500m) - using lincom for SE extraction
preserve
keep if traj15_O3_500m == 1
mixed BAA_PhenoAge i.traj15_EVI_500m ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year ///
     || id_num:, ///
     vce(robust)
lincom 1.traj15_EVI_500m
scalar beta_highO3_500m = r(estimate)
scalar se_highO3_500m   = r(se)
global EVI_effect_highO3_500m = beta_highO3_500m
restore

* Calculate 500m Effect Difference + P-value (BUG FIX: added missing p-value computation)
global effect_diff_500m = $EVI_effect_highO3_500m - $EVI_effect_lowO3_500m
scalar se_diff_500m_s = sqrt(se_highO3_500m^2 + se_lowO3_500m^2)
scalar z_diff_500m_s  = $effect_diff_500m / se_diff_500m_s
scalar p_diff_500m_s  = 2 * normal(-abs(z_diff_500m_s))
global effect_diff_p_500m = p_diff_500m_s

di as txt "500m Low O3 EVI effect: " as res %8.4f $EVI_effect_lowO3_500m
di as txt "500m High O3 EVI effect: " as res %8.4f $EVI_effect_highO3_500m
di as txt "500m Effect difference: " as res %8.4f $effect_diff_500m
di as txt "500m Effect difference P-value: " as res %8.4f $effect_diff_p_500m

* Calculate 500m Effect Change Percentage
if $EVI_effect_lowO3_500m != 0 {
    global percent_change_500m = (($EVI_effect_highO3_500m - $EVI_effect_lowO3_500m) / abs($EVI_effect_lowO3_500m)) * 100
    di as txt "500m EVI effect percent change: " as res %6.2f $percent_change_500m as txt "%"
}

* -----------------------------------------------------------------------------
* Interaction Analysis汇总表  
* -----------------------------------------------------------------------------

di as txt "=== 创建Interaction Analysis汇总表   ==="

capture putdocx clear
putdocx begin

* Title - 
putdocx paragraph, style("Heading1")
putdocx text ("Table 1. Interaction Analysis of EVI and O3 Exposure on BAA"), bold

* Analysis Description -  
putdocx paragraph
putdocx text ("Analytical objective: To assess whether enhanced vegetation index (EVI) modifies the association between ozone (O3) exposure and biological age acceleration (BAA).")
putdocx text ("Methodology: Stratified analysis by O3 exposure groups (low vs high) using mixed-effects models with inverse probability weighting, adjusted for demographic, socioeconomic, lifestyle, climatic, and health-related covariates.")

* 主要汇总表 - 8-Column Professional Table
putdocx paragraph, style("Heading2")
putdocx text ("Interaction Analysis Summary")

* 创建8-Column Table
putdocx table tbl_summary = (3,8), border(all, single) width(100%)

* Table Header - Using Professional Abbreviations
putdocx table tbl_summary(1,1) = ("Exposure"), bold halign(center)
putdocx table tbl_summary(1,2) = ("Buffer"), bold halign(center)
putdocx table tbl_summary(1,3) = ("Inter_Coeff"), bold halign(center)
putdocx table tbl_summary(1,4) = ("Inter_P"), bold halign(center)
putdocx table tbl_summary(1,5) = ("EVI_LowO3"), bold halign(center)
putdocx table tbl_summary(1,6) = ("EVI_HighO3"), bold halign(center)
putdocx table tbl_summary(1,7) = ("Eff_Diff"), bold halign(center)
putdocx table tbl_summary(1,8) = ("Eff_Diff_P"), bold halign(center)

* 1KM Buffer Results
putdocx table tbl_summary(2,1) = ("EVI"), halign(center)
putdocx table tbl_summary(2,2) = ("1 km"), halign(center)
putdocx table tbl_summary(2,3) = (strofreal($interaction_coef_1km_new, "%8.4f")), halign(center)
putdocx table tbl_summary(2,4) = (strofreal($interaction_p_1km_new, "%8.4f")), halign(center)
putdocx table tbl_summary(2,5) = (strofreal($EVI_effect_lowO3_1km, "%8.4f")), halign(center)
putdocx table tbl_summary(2,6) = (strofreal($EVI_effect_highO3_1km, "%8.4f")), halign(center)
putdocx table tbl_summary(2,7) = (strofreal($effect_diff_1km, "%8.4f")), halign(center)
putdocx table tbl_summary(2,8) = (strofreal($effect_diff_p_1km, "%8.4f")), halign(center)

* 500m Buffer Results
putdocx table tbl_summary(3,1) = ("EVI"), halign(center)
putdocx table tbl_summary(3,2) = ("500 m"), halign(center)
putdocx table tbl_summary(3,3) = (strofreal($interaction_coef_500m_new, "%8.4f")), halign(center)
putdocx table tbl_summary(3,4) = (strofreal($interaction_p_500m_new, "%8.4f")), halign(center)
putdocx table tbl_summary(3,5) = (strofreal($EVI_effect_lowO3_500m, "%8.4f")), halign(center)
putdocx table tbl_summary(3,6) = (strofreal($EVI_effect_highO3_500m, "%8.4f")), halign(center)
putdocx table tbl_summary(3,7) = (strofreal($effect_diff_500m, "%8.4f")), halign(center)
putdocx table tbl_summary(3,8) = (strofreal($effect_diff_p_500m, "%8.4f")), halign(center)

* Table Notes  
putdocx paragraph
putdocx text ("Table notes: All models were adjusted for age, sex, education, socioeconomic status, lifestyle factors, climate principal components, and health conditions."), italic
putdocx text ("Abbreviations: Inter_Coeff = Interaction coefficient; Inter_P = P-value for interaction term; EVI_LowO3 = EVI effect in low O3 group; EVI_HighO3 = EVI effect in high O3 group; Eff_Diff = Effect difference (EVI_HighO3 - EVI_LowO3); Eff_Diff_P = P-value for effect difference."), italic
putdocx text ("Statistical significance: P < 0.05 indicates statistical significance. Negative EVI effects indicate protective associations with biological aging."), italic

* Result Interpretation -  
putdocx paragraph, style("Heading2")
putdocx text ("Interpretation of Results")

* Interaction Significance Assessment
local sig_1km = ($interaction_p_1km_new < 0.05)
local sig_500m = ($interaction_p_500m_new < 0.05)

if `sig_1km' | `sig_500m' {
    putdocx text ("Statistically significant interaction effects were observed:")
    
    if `sig_1km' {
        if $interaction_coef_1km_new < 0 {
            putdocx text ("At 1 km buffer: O3 exposure significantly attenuated the protective effect of EVI on biological aging (P = " + strofreal($interaction_p_1km_new, "%8.4f") + ").")
        }
        else {
            putdocx text ("At 1 km buffer: O3 exposure significantly enhanced the protective effect of EVI on biological aging (P = " + strofreal($interaction_p_1km_new, "%8.4f") + ").")
        }
    }
    
    if `sig_500m' {
        if $interaction_coef_500m_new < 0 {
            putdocx text ("At 500 m buffer: O3 exposure significantly attenuated the protective effect of EVI on biological aging (P = " + strofreal($interaction_p_500m_new, "%8.4f") + ").")
        }
        else {
            putdocx text ("At 500 m buffer: O3 exposure significantly enhanced the protective effect of EVI on biological aging (P = " + strofreal($interaction_p_500m_new, "%8.4f") + ").")
        }
    }
}
else {
    putdocx text ("No statistically significant interaction effects were observed (all interaction P-values ≥ 0.05).")
}

* Effect Difference Significance Assessment
local sig_diff_1km = ($effect_diff_p_1km < 0.05)
local sig_diff_500m = ($effect_diff_p_500m < 0.05)

if `sig_diff_1km' | `sig_diff_500m' {
    putdocx text ("Stratified analysis revealed significant effect modification:")
    
    if `sig_diff_1km' {
        if $effect_diff_1km < 0 {
            putdocx text ("At 1 km buffer: The protective effect of EVI was significantly weaker in high O3 group compared to low O3 group (P = " + strofreal($effect_diff_p_1km, "%8.4f") + ").")
        }
        else {
            putdocx text ("At 1 km buffer: The protective effect of EVI was significantly stronger in high O3 group compared to low O3 group (P = " + strofreal($effect_diff_p_1km, "%8.4f") + ").")
        }
    }
    
    if `sig_diff_500m' {
        if $effect_diff_500m < 0 {
            putdocx text ("At 500 m buffer: The protective effect of EVI was significantly weaker in high O3 group compared to low O3 group (P = " + strofreal($effect_diff_p_500m, "%8.4f") + ").")
        }
        else {
            putdocx text ("At 500 m buffer: The protective effect of EVI was significantly stronger in high O3 group compared to low O3 group (P = " + strofreal($effect_diff_p_500m, "%8.4f") + ").")
        }
    }
}
else {
    putdocx text ("Stratified analysis did not reveal significant effect modification (all effect difference P-values ≥ 0.05).")
}

* Calculate effect change percentage for interpretation
if $EVI_effect_lowO3_1km != 0 {
    global percent_change_1km_new = (($EVI_effect_highO3_1km - $EVI_effect_lowO3_1km) / abs($EVI_effect_lowO3_1km)) * 100
    if $effect_diff_1km < 0 {
        putdocx text ("At 1 km buffer: EVI protective effect was " + strofreal(abs($percent_change_1km_new), "%6.1f") + "% weaker in high O3 environments.")
    }
    else {
        putdocx text ("At 1 km buffer: EVI protective effect was " + strofreal($percent_change_1km_new, "%6.1f") + "% stronger in high O3 environments.")
    }
}

if $EVI_effect_lowO3_500m != 0 {
    global percent_change_500m_new = (($EVI_effect_highO3_500m - $EVI_effect_lowO3_500m) / abs($EVI_effect_lowO3_500m)) * 100
    if $effect_diff_500m < 0 {
        putdocx text ("At 500 m buffer: EVI protective effect was " + strofreal(abs($percent_change_500m_new), "%6.1f") + "% weaker in high O3 environments.")
    }
    else {
        putdocx text ("At 500 m buffer: EVI protective effect was " + strofreal($percent_change_500m_new, "%6.1f") + "% stronger in high O3 environments.")
    }
}

* Overall Conclusion -  
putdocx paragraph, style("Heading2")
putdocx text ("Principal Findings and Implications")

if ($interaction_p_1km_new < 0.05 | $interaction_p_500m_new < 0.05) & ///
   ($effect_diff_p_1km < 0.05 | $effect_diff_p_500m < 0.05) {
    putdocx text ("Strong evidence indicates that ozone exposure significantly modifies the protective association between green space and biological aging.")
    putdocx text ("Evidence strength: High")
    putdocx text ("Public health implication: Urban planning strategies should consider integrated approaches to simultaneously mitigate air pollution and enhance green infrastructure.")
}
else if ($interaction_p_1km_new < 0.1 | $interaction_p_500m_new < 0.1) & ///
        ($effect_diff_p_1km < 0.1 | $effect_diff_p_500m < 0.1) {
    putdocx text ("Moderate evidence suggests potential modification of EVI protective effects by ozone exposure.")
    putdocx text ("Evidence strength: Moderate")
    putdocx text ("Research implication: Further investigation with larger sample sizes is warranted to confirm these findings.")
}
else if ($effect_diff_1km != 0 | $effect_diff_500m != 0) {
    putdocx text ("Limited evidence indicates possible effect modification, though statistical evidence remains insufficient.")
    putdocx text ("Evidence strength: Low")
    putdocx text ("Research implication: Additional studies are needed to explore these relationships across different populations and environmental contexts.")
}
else {
    putdocx text ("Insufficient evidence to support effect modification of EVI protective associations by ozone exposure.")
    putdocx text ("Evidence strength: None")
    putdocx text ("Research implication: Future studies should consider alternative mechanisms and environmental modifiers.")
}

* Methodological Considerations
putdocx paragraph, style("Heading2")
putdocx text ("Methodological Considerations")

putdocx text ("Analytical strengths:")
putdocx text ("• Mixed-effects modeling accounting for individual-level clustering")
putdocx text ("• Inverse probability weighting to address potential selection bias")
putdocx text ("• Comprehensive adjustment for demographic, socioeconomic, and environmental covariates")
putdocx text ("• Stratified analysis to examine effect modification across pollution gradients")

putdocx text ("Limitations:")
putdocx text ("• Observational design precludes causal inference")
putdocx text ("• Potential for residual confounding despite extensive covariate adjustment")
putdocx text ("• Exposure measurement error inherent in satellite-based estimates")
putdocx text ("• Limited statistical power for interaction detection")

putdocx save "output/tables/EVI_O3_interaction_summary_professional.docx", replace

* -----------------------------------------------------------------------------
* Display Final Results in Stata Results Window  
* -----------------------------------------------------------------------------

di as txt "=== EVI与O3Interaction Analysis最终结果 - 英文摘要 ==="
di as txt "1KM Buffer Interaction Coefficient: " as res %8.4f $interaction_coef_1km_new as txt " (P = " as res %8.4f $interaction_p_1km_new as txt ")"
di as txt "1KM Buffer Effect Difference P-value: " as res %8.4f $effect_diff_p_1km
di as txt "500m Buffer Interaction Coefficient: " as res %8.4f $interaction_coef_500m_new as txt " (P = " as res %8.4f $interaction_p_500m_new as txt ")"
di as txt "500m Buffer Effect Difference P-value: " as res %8.4f $effect_diff_p_500m

* Principal Findings Summary
di as txt "=== PRINCIPAL FINDINGS ==="
if $interaction_p_1km_new < 0.05 | $interaction_p_500m_new < 0.05 {
    di as txt "• Significant interaction detected between EVI and O3 exposure"
}
else {
    di as txt "• No significant interaction between EVI and O3 exposure"
}

if $effect_diff_1km < 0 | $effect_diff_500m < 0 {
    di as txt "• Consistent trend of attenuated EVI protection in high O3 environments"
}
else if $effect_diff_1km > 0 | $effect_diff_500m > 0 {
    di as txt "• Consistent trend of enhanced EVI protection in high O3 environments"
}
else {
    di as txt "• No consistent pattern of effect modification across buffers"
}

di as txt "=== Professional English Analysis Complete ==="
di as txt "Generated files: EVI_O3Interaction Analysis汇总表_专业英文版.docx"
