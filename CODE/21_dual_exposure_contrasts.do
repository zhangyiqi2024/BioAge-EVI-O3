* =============================================================================
* Analysis of Partial Offsetting of O3 Damage Effect by EVI Based on Joint Effect Model Results
* =============================================================================

set seed 20240101

di as txt "=== Analysis of Partial Offsetting of O3 Damage by EVI Based on Joint Effect Models ==="

* Reload data
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

* Principal component analysis of climate variables
di "=== Principal Component Analysis of Climate Variables ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* -----------------------------------------------------------------------------
* Re-run four joint effect models and extract key results
* -----------------------------------------------------------------------------

* Model 3-1: EVI_1KM + O3_1KM Joint Effect on BAA_PhenoAge (Mixed-effects Model)
di as txt "=== Model 3-1: Joint Effect Analysis at 1KM Scale ==="
mixed BAA_PhenoAge i.joint_exp_1KM                   ///
     c.age i.sex i.education                         ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current              ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia              ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_1KM) vce(robust)

estimates store joint_1km
matrix b_1km = e(b)
matrix V_1km = e(V)

* Extract key group comparisons
di as txt "Key Group Comparisons (1KM Scale):"

* 1. O3 individual damage effect: Low EVI High O3 vs Low EVI Low O3
lincom 2.joint_exp_1KM - 1.joint_exp_1KM
local O3_effect_1km = r(estimate)
local O3_p_1km = r(p)
di as txt "  O3 Individual Damage Effect = " as res %8.4f `O3_effect_1km' as txt " (p = " as res %8.4f `O3_p_1km' as txt ")"

* 2. EVI individual protective effect: High EVI Low O3 vs Low EVI Low O3
lincom 3.joint_exp_1KM - 1.joint_exp_1KM
local EVI_effect_1km = r(estimate)
local EVI_p_1km = r(p)
di as txt "  EVI Individual Protective Effect = " as res %8.4f `EVI_effect_1km' as txt " (p = " as res %8.4f `EVI_p_1km' as txt ")"

* 3. Joint exposure net effect: High EVI High O3 vs Low EVI Low O3
lincom 4.joint_exp_1KM - 1.joint_exp_1KM
local joint_effect_1km = r(estimate)
local joint_p_1km = r(p)
di as txt "  Joint Exposure Net Effect = " as res %8.4f `joint_effect_1km' as txt " (p = " as res %8.4f `joint_p_1km' as txt ")"

* 4. Offset effect: High EVI High O3 vs Low EVI High O3
lincom 4.joint_exp_1KM - 2.joint_exp_1KM
local offset_effect_1km = r(estimate)
local offset_p_1km = r(p)
di as txt "  Offset Effect = " as res %8.4f `offset_effect_1km' as txt " (p = " as res %8.4f `offset_p_1km' as txt ")"

* Model 3-2: EVI_500m + O3_500m Joint Effect on BAA_PhenoAge (Mixed-effects Model)
di as txt "=== Model 3-2: Joint Effect Analysis at 500m Scale ==="
mixed BAA_PhenoAge i.joint_exp_500m                  ///
     c.age i.sex i.education                         ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current              ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia              ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_500m) vce(robust)

estimates store joint_500m
matrix b_500m = e(b)
matrix V_500m = e(V)

* Extract key group comparisons
di as txt "Key Group Comparisons (500m Scale):"

lincom 2.joint_exp_500m - 1.joint_exp_500m
local O3_effect_500m = r(estimate)
local O3_p_500m = r(p)

lincom 3.joint_exp_500m - 1.joint_exp_500m
local EVI_effect_500m = r(estimate)
local EVI_p_500m = r(p)

lincom 4.joint_exp_500m - 1.joint_exp_500m
local joint_effect_500m = r(estimate)
local joint_p_500m = r(p)

lincom 4.joint_exp_500m - 2.joint_exp_500m
local offset_effect_500m = r(estimate)
local offset_p_500m = r(p)

di as txt "  O3 Individual Damage Effect = " as res %8.4f `O3_effect_500m' as txt " (p = " as res %8.4f `O3_p_500m' as txt ")"
di as txt "  EVI Individual Protective Effect = " as res %8.4f `EVI_effect_500m' as txt " (p = " as res %8.4f `EVI_p_500m' as txt ")"
di as txt "  Joint Exposure Net Effect = " as res %8.4f `joint_effect_500m' as txt " (p = " as res %8.4f `joint_p_500m' as txt ")"
di as txt "  Offset Effect = " as res %8.4f `offset_effect_500m' as txt " (p = " as res %8.4f `offset_p_500m' as txt ")"

* -----------------------------------------------------------------------------
* Calculate offset magnitude and evidence for partial offsetting
* -----------------------------------------------------------------------------
di as txt "=== Evidence Analysis for Partial Offsetting of O3 Damage by EVI ==="

* 1. Calculate theoretical offset percentage
if `O3_effect_1km' != 0 {
    local offset_percentage_1km = (`offset_effect_1km' / `O3_effect_1km') * 100
    di as txt "Theoretical Offset Percentage at 1KM Scale = " as res %8.1f `offset_percentage_1km' as txt "%"
}

if `O3_effect_500m' != 0 {
    local offset_percentage_500m = (`offset_effect_500m' / `O3_effect_500m') * 100
    di as txt "Theoretical Offset Percentage at 500m Scale = " as res %8.1f `offset_percentage_500m' as txt "%"
}

* 2. Assess evidence for partial offsetting
di as txt "Evidence Criteria for Partial Offsetting:"
di as txt "1. EVI individual effect is negative (protective)"
di as txt "2. O3 individual effect is positive (damaging)"
di as txt "3. Joint exposure net effect lies between EVI and O3 individual effects"
di as txt "4. Offset effect direction is consistent with EVI protection"

local evidence_count = 0

* Check evidence at 1KM scale
di as txt "Evidence at 1KM Scale:"
if `EVI_effect_1km' < 0 {
    di as txt "  ✓ EVI individual effect is negative (protective)"
    local evidence_count = `evidence_count' + 1
}

if `O3_effect_1km' > 0 {
    di as txt "  ✓ O3 individual effect is positive (damaging)"
    local evidence_count = `evidence_count' + 1
}

if `joint_effect_1km' > `EVI_effect_1km' & `joint_effect_1km' < `O3_effect_1km' {
    di as txt "  ✓ Joint exposure net effect lies between EVI protection and O3 damage effects"
    local evidence_count = `evidence_count' + 1
}

if `offset_effect_1km' < 0 {
    di as txt "  ✓ Offset effect direction is consistent with EVI protection"
    local evidence_count = `evidence_count' + 1
}

* Check evidence at 500m scale
local evidence_count_500m = 0
di as txt "Evidence at 500m Scale:"
if `EVI_effect_500m' < 0 {
    di as txt "  ✓ EVI individual effect is negative (protective)"
    local evidence_count_500m = `evidence_count_500m' + 1
}

if `O3_effect_500m' > 0 {
    di as txt "  ✓ O3 individual effect is positive (damaging)"
    local evidence_count_500m = `evidence_count_500m' + 1
}

if `joint_effect_500m' > `EVI_effect_500m' & `joint_effect_500m' < `O3_effect_500m' {
    di as txt "  ✓ Joint exposure net effect lies between EVI protection and O3 damage effects"
    local evidence_count_500m = `evidence_count_500m' + 1
}

if `offset_effect_500m' < 0 {
    di as txt "  ✓ Offset effect direction is consistent with EVI protection"
    local evidence_count_500m = `evidence_count_500m' + 1
}

* -----------------------------------------------------------------------------
* Generate visualization evidence plot
* -----------------------------------------------------------------------------
di as txt "=== Generating Joint Exposure Group Effect Pattern Plot ==="

* Create plot data
preserve
clear
set obs 8

generate scale = ""
generate group = ""
generate effect = .
generate group_num = .

* 1KM scale data
replace scale = "1KM" if _n <= 4
replace group_num = _n if _n <= 4
replace group = "Low EVI Low O3" if _n == 1
replace group = "Low EVI High O3" if _n == 2
replace group = "High EVI Low O3" if _n == 3
replace group = "High EVI High O3" if _n == 4

replace effect = 0 if _n == 1  // Reference group
replace effect = `O3_effect_1km' if _n == 2
replace effect = `EVI_effect_1km' if _n == 3
replace effect = `joint_effect_1km' if _n == 4

* 500m scale data
replace scale = "500m" if _n > 4
replace group_num = _n - 4 if _n > 4
replace group = "Low EVI Low O3" if _n == 5
replace group = "Low EVI High O3" if _n == 6
replace group = "High EVI Low O3" if _n == 7
replace group = "High EVI High O3" if _n == 8

replace effect = 0 if _n == 5  // Reference group
replace effect = `O3_effect_500m' if _n == 6
replace effect = `EVI_effect_500m' if _n == 7
replace effect = `joint_effect_500m' if _n == 8

* Generate plot position variable
generate x_pos = group_num + (scale == "500m") * 0.4 - 0.2

* Draw effect pattern plot (name it for later combining)
twoway (bar effect x_pos if scale == "1KM", barw(0.2) color(blue%70)) ///
       (bar effect x_pos if scale == "500m", barw(0.2) color(red%70)) ///
       (scatteri 0.2 1 "Reference", msymbol(none) mlabel(Reference) mlabpos(12) mlabsize(small) mlabcolor(black)), ///
       xlabel(1 `"Low EVI Low O{subscript:3}"' 2 `"Low EVI High O{subscript:3}"' 3 `"High EVI Low O{subscript:3}"' 4 `"High EVI High O{subscript:3}"', labsize(small)) ///
       ytitle("BAA difference vs. Low EVI & Low O3 (years)", size(small)) ///
       xtitle("") ///
       legend(order(1 "1km scale" 2 "500m scale") rows(1) size(small)) ///
       xscale(range(0.5, 4.5)) ///
       yline(0, lcolor(black) lpattern(dash)) ///
       plotregion(margin(zero)) ///
       title("b", size(large) position(11) ring(0)) ///
       note(`"Note: Positive values indicate biological age acceleration (harmful), negative values indicate biological age deceleration (beneficial)"' ///
            `"The effect of High EVI High O{subscript:3} group lies between High EVI Low O{subscript:3} and Low EVI High O{subscript:3} groups, suggesting apparent additive offsetting."', size(vsmall)) ///
       name(offset_plot, replace)

graph export "output/figures/offset_effect_pattern.png", replace width(2000) height(1200)

* ---- Save offset graph to disk ----
graph save offset_plot "output/figures/_offset_plot.gph", replace
restore

* ---- Load joint exposure graph (from 07_joint_exposure_groups.do) and combine ----
graph use "output/figures/_joint_exp.gph", name(joint_exp, replace)
graph use "output/figures/_offset_plot.gph", name(offset_plot, replace)

graph combine joint_exp offset_plot, ///
    rows(2) title("") ysize(10) xsize(8) ///
    graphregion(color(white))

graph export "output/figures/Figure3_combined_joint_offset.png", replace width(2000) height(2400)
erase "output/figures/_joint_exp.gph"
erase "output/figures/_offset_plot.gph"
di as txt "✓ Combined Figure 3 saved: ../output/figures/Figure3_combined_joint_offset.png"

* -----------------------------------------------------------------------------
* Comprehensive conclusion
* -----------------------------------------------------------------------------
di as txt "=== Comprehensive Conclusion on Partial Offsetting of O3 Damage by EVI ==="

di as txt "Evidence Summary:"
di as txt "1KM Scale Evidence Criteria Met: " as res `evidence_count' as txt "/4"
di as txt "500m Scale Evidence Criteria Met: " as res `evidence_count_500m' as txt "/4"

* Draw conclusion based on evidence
if `evidence_count' >= 3 | `evidence_count_500m' >= 3 {
    di as txt "Conclusion: " as res "✅ Clear evidence exists for partial offsetting of O3 damage by EVI"
    di as txt "      - EVI protective effect partially offsets O3 damage effect"
    di as txt "      - High vegetation coverage environments can mitigate the harm of high ozone exposure"
}
else if `evidence_count' >= 2 | `evidence_count_500m' >= 2 {
    di as txt "Conclusion: " as res "⚠️ Suggestive evidence exists for partial offsetting of O3 damage by EVI"
    di as txt "      - Observed pattern suggests EVI may partially mitigate O3 harm"
    di as txt "      - More evidence needed to confirm this relationship"
}
else {
    di as txt "Conclusion: " as res "❓ Current evidence insufficient to support partial offsetting of O3 damage by EVI"
    di as txt "      - No consistent offsetting pattern observed across joint exposure groups"
}

* Generate detailed report
capture putdocx clear
putdocx begin
putdocx paragraph, style("Heading1")
putdocx text ("Analysis Report on Partial Offsetting of O3 Damage by EVI")

putdocx paragraph
putdocx text ("Based on results from four joint exposure models, systematically evaluate whether vegetation coverage (EVI) partially offsets the damage effect of ozone (O3) exposure on biological age acceleration.")

putdocx paragraph, style("Heading2")
putdocx text ("Key Findings")

putdocx paragraph
putdocx text ("1. EVI Individual Effect (Protective):")
putdocx text ("   - 1KM Scale: ")
putdocx text ("`EVI_effect_1km'")
putdocx text (" (p=")
putdocx text ("`EVI_p_1km'")
putdocx text (")")
putdocx text ("   - 500m Scale: ")
putdocx text ("`EVI_effect_500m'")
putdocx text (" (p=")
putdocx text ("`EVI_p_500m'")
putdocx text (")")

putdocx paragraph
putdocx text ("2. O3 Individual Effect (Damaging):")
putdocx text ("   - 1KM Scale: ")
putdocx text ("`O3_effect_1km'")
putdocx text (" (p=")
putdocx text ("`O3_p_1km'")
putdocx text (")")
putdocx text ("   - 500m Scale: ")
putdocx text ("`O3_effect_500m'")
putdocx text (" (p=")
putdocx text ("`O3_p_500m'")
putdocx text (")")

putdocx paragraph
putdocx text ("3. Joint Exposure Net Effect:")
putdocx text ("   - 1KM Scale: ")
putdocx text ("`joint_effect_1km'")
putdocx text (" (p=")
putdocx text ("`joint_p_1km'")
putdocx text (")")
putdocx text ("   - 500m Scale: ")
putdocx text ("`joint_effect_500m'")
putdocx text (" (p=")
putdocx text ("`joint_p_500m'")
putdocx text (")")

putdocx paragraph
putdocx text ("4. Offset Effect:")
putdocx text ("   - 1KM Scale: ")
putdocx text ("`offset_effect_1km'")
putdocx text (" (p=")
putdocx text ("`offset_p_1km'")
putdocx text (")")
putdocx text ("   - 500m Scale: ")
putdocx text ("`offset_effect_500m'")
putdocx text (" (p=")
putdocx text ("`offset_p_500m'")
putdocx text (")")

if `O3_effect_1km' != 0 {
    putdocx paragraph
    putdocx text ("5. Theoretical Offset Percentage:")
    putdocx text ("   - 1KM Scale: ")
    putdocx text ("`offset_percentage_1km'%")
    if `O3_effect_500m' != 0 {
        putdocx text ("   - 500m Scale: ")
        putdocx text ("`offset_percentage_500m'%")
    }
}

putdocx paragraph, style("Heading2")
putdocx text ("Comprehensive Conclusion")
if `evidence_count' >= 3 | `evidence_count_500m' >= 3 {
    putdocx text ("Based on joint exposure model results, clear evidence exists that EVI partially offsets the damage effect of O3 on biological age acceleration.")
    putdocx text ("High vegetation coverage environments can mitigate the health hazards of high ozone exposure to some extent.")
}
else if `evidence_count' >= 2 | `evidence_count_500m' >= 2 {
    putdocx text ("Based on joint exposure model results, suggestive evidence exists that EVI may partially offset the damage effect of O3.")
    putdocx text ("Further research is needed to confirm this relationship.")
}
else {
    putdocx text ("Based on current joint exposure model results, evidence is insufficient to support the conclusion that EVI partially offsets O3 damage effect.")
}

putdocx save "output/tables/offset_analysis_report.docx", replace

di as txt "=== Analysis Completed ==="
di as txt "Analysis report saved as: ../output/tables/offset_analysis_report.docx"
di as txt "Effect pattern plot saved as: ../output/figures/offset_effect_pattern.png"
