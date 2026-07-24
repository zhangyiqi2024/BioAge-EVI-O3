*------------------------------------------------------------
* Dual-exposure model counterfactual prediction analysis: Main model + counterfactual prediction + by age group + figure export
*------------------------------------------------------------

set seed 20240101

* Create output folders
capture mkdir "output/figures/counterfactual/"
capture mkdir "output/tables/"

* Load data
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90

gen ln_hh_consumption_total = log(hh_consumption_total)

* Climate variable principal component analysis
di "=== Climate PCA ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* Run dual-exposure main effects model
mixed BAA_PhenoAge c.EVI_1KM_mean4y c.O3_1KM_mean4y ///
     c.age i.sex i.education ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance ///
     i.smoking_current i.alcohol_current ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
     c.bmi i.hypertension i.dyslipidemia ///
     i.year ///
     || id:, ///
     vce(robust)

* Save model results
estimates store double_exposure_model_main

* Display model results
estimates replay double_exposure_model_main

* Extract model coefficients
matrix beta = e(b)
local evi_coef = beta[1,1]
local o3_coef  = beta[1,2]

display "EVI coefficient: `evi_coef'"
display "O3 coefficient: `o3_coef'"

* Predict actual BAA_PhenoAge (fixed effects component)
predict actual_BAA_fixed, xb

* Create counterfactual exposure variables
gen counterfactual_EVI = EVI_1KM_mean4y
replace counterfactual_EVI = 0.40 if EVI_1KM_mean4y < 0.40

gen counterfactual_O3  = O3_1KM_mean4y
replace counterfactual_O3 = 60 if O3_1KM_mean4y > 60

* Manually calculate counterfactual predicted values (fixed effects component)
gen counterfactual_BAA_fixed = actual_BAA_fixed - ///
     (`evi_coef' * EVI_1KM_mean4y) - (`o3_coef' * O3_1KM_mean4y) + ///
     (`evi_coef' * counterfactual_EVI) + (`o3_coef' * counterfactual_O3)

* Calculate BAA_PhenoAge reduction (positive values indicate improvement)
gen BAA_reduction = actual_BAA_fixed - counterfactual_BAA_fixed

* Analyze exposure improvement
gen EVI_improvement = counterfactual_EVI - EVI_1KM_mean4y
gen O3_improvement  = counterfactual_O3  - O3_1KM_mean4y

* Define improvement type groups
gen improvement_type = ""
replace improvement_type = "both_improved"      if EVI_improvement > 0 & O3_improvement < 0
replace improvement_type = "only_evi_improved"  if EVI_improvement > 0 & O3_improvement >= 0
replace improvement_type = "only_o3_improved"   if EVI_improvement <= 0 & O3_improvement < 0
replace improvement_type = "no_improvement"     if EVI_improvement <= 0 & O3_improvement >= 0

* Create age group variable (in correct order)
gen age_group = ""
replace age_group = "<=60 years"   if age <= 60
replace age_group = "61-70 years" if age > 60 & age <= 70
replace age_group = "71-80 years" if age > 70 & age <= 80
replace age_group = ">80 years"   if age > 80 & age < 90

* Create age group order variable to ensure correct ordering
gen age_group_order = .
replace age_group_order = 1 if age_group == "<=60 years"
replace age_group_order = 2 if age_group == "61-70 years"
replace age_group_order = 3 if age_group == "71-80 years"
replace age_group_order = 4 if age_group == ">80 years"

label define age_group_label 1 "<=60 years" 2 "61-70 years" 3 "71-80 years" 4 ">80 years"
label values age_group_order age_group_label

* Check group distribution
display "=== Improvement type distribution ==="
tabulate improvement_type

* Basic descriptive statistics
display "=== Overall description of BAA_PhenoAge reduction ==="
summarize BAA_reduction, detail

display "=== Mean and 95% CI of BAA_PhenoAge reduction ==="
ci mean BAA_reduction

* Calculate improvement proportion
gen baa_improved = BAA_reduction > 0
label variable baa_improved "BAA improved (yes=1)"
tabstat baa_improved, stat(mean) format(%9.3f)
display "Improvement proportion: " %4.1f (r(mean)*100) "%"

* Create valid data subset for subsequent analysis
tempfile valid_data
save "`valid_data'", replace

use "`valid_data'", clear
keep if !missing(BAA_reduction)

* Check improvement type distribution in valid data
display "=== Improvement type distribution in valid data ==="
tabulate improvement_type

* BAA reduction by improvement type
display "=== BAA reduction by improvement type ==="
tabstat BAA_reduction, by(improvement_type) stat(mean sd n)

* Check group sample sizes
count if improvement_type == "both_improved"
local both_improved_n = r(N)

count if improvement_type == "only_evi_improved"
local only_evi_improved_n = r(N)

count if improvement_type == "only_o3_improved"
local only_o3_improved_n = r(N)

count if improvement_type == "no_improvement"
local no_improvement_n = r(N)

display "Valid sample size by group:"
display "Both improved: `both_improved_n'"
display "Only EVI improved: `only_evi_improved_n'"
display "Only O3 improved: `only_o3_improved_n'"
display "No improvement: `no_improvement_n'"

* Between-group difference tests - only for valid groups
if `both_improved_n' > 0 & `only_o3_improved_n' > 0 {
    display "=== Between-group difference tests ==="

    * Both improved vs Only O3 improved
    ttest BAA_reduction if improvement_type == "both_improved" | improvement_type == "only_o3_improved", by(improvement_type)
    display "Both improved vs Only O3 improved:"
    display "  t = " %6.3f r(t) ", P = " %6.4f r(p)
    display "  Both improved mean = " %6.3f r(mu_1) " years"
    display "  Only O3 improved mean = " %6.3f r(mu_2) " years"
    display "  Mean difference = " %6.3f (r(mu_1) - r(mu_2)) " years"
}

* Check age group distribution
display "=== Age group distribution ==="
tabulate age_group

*------------------------------------------------------------
* Figure 1: BAA improvement box plot by age group
*------------------------------------------------------------
graph box BAA_reduction, over(age_group_order, label(labsize(small))) ///
    title("BAA improvement by age group", size(medium)) ///
    ytitle("BAA improvement (years)", size(small)) ///
    marker(1, msymbol(oh)) ///
    box(1, color("30 120 180")) ///
    intensity(70) ///
    note("Positive values indicate BAA improvement (lower biological age)", size(vsmall))
graph export "output/figures/counterfactual/baa_improvement_boxplot_by_age.png", replace width(2000) height(1500)

******************************************************
* Add "All participants" to the BAA improvement box plot by age group and improvement type
* Use English throughout for figures and CSVs, and reduce x-axis label size to avoid overlap
******************************************************
preserve
    * Use valid data
    use "`valid_data'", clear
    keep if !missing(BAA_reduction)

    * Extended age group variable (including all participants), English labels
    gen age_group_ext = age_group_order
    capture label drop agegroupl
    label define agegroupl ///
        1 "<=60 years" ///
        2 "61-70 years" ///
        3 "71-80 years" ///
        4 ">80 years" ///
        5 "All participants"
    label values age_group_ext agegroupl

    * English version of improvement type (for plotting and export)
    gen improve_type_en = ""
    replace improve_type_en = "Both improved" if improvement_type == "both_improved"
    replace improve_type_en = "Only EVI improved"                    if improvement_type == "only_evi_improved"
    replace improve_type_en = "Only O3 improved"                     if improvement_type == "only_o3_improved"
    replace improve_type_en = "No improvement"                       if improvement_type == "no_improvement"

    * Save age-specific raw data
    tempfile age_specific overall_data
    save `age_specific'

    * Generate "All participants (45<age<90)" data, keeping only "Only O3 improved" and "Both improved"
    keep if inlist(improvement_type, "only_o3_improved", "both_improved")
    replace age_group_ext = 5
    save `overall_data', replace

    * Restore age-stratified data and append all-participants
    use `age_specific', clear
    append using `overall_data'

    * Re-generate English improvement type variable to ensure appended portion also has it
    drop improve_type_en
    gen improve_type_en = ""
    replace improve_type_en = "Both improved" if improvement_type == "both_improved"
    replace improve_type_en = "Only EVI improved"                    if improvement_type == "only_evi_improved"
    replace improve_type_en = "Only O3 improved"                     if improvement_type == "only_o3_improved"
    replace improve_type_en = "No improvement"                       if improvement_type == "no_improvement"

    * Draw English box plot:
    *   Outer: Age group (age_group_ext)
    *   Inner: Improvement type (improve_type_en)
    *   X-axis improvement type labels: vsmall + angle(45), to avoid overlap
    graph box BAA_reduction, ///
        over(improve_type_en, label(labsize(vsmall) angle(45))) ///
        over(age_group_ext, label(labsize(small))) ///
        title("BAA improvement by age group and exposure improvement pattern", ///
              size(medium)) ///
        ytitle("BAA improvement (years)", size(small)) ///
        marker(1, msymbol(oh)) ///
        box(1, color("30 120 180")) ///
        intensity(70) ///
        note("Positive values indicate improvement (lower biological age). 'All participants' refers to 45<age<90.", ///
             size(vsmall))

    * Overwrite original filename (figure content is entirely in English)
    graph export "output/figures/counterfactual/baa_improvement_boxplot_by_age_and_type.png", ///
        replace width(2000) height(1500)

    **************************************************
    * Export statistical parameters for this figure to CSV (all English)
    **************************************************
    levelsof age_group_ext,   local(age_levels)
    levelsof improve_type_en, local(type_levels)

    tempname memhold
    tempfile box_params

    * Coeff = mean, StdErr = SE(mean), P_value = H0: mean = 0, CI_* = 95%CI
    postfile `memhold' ///
        str40 AgeGroup str50 ImprovementType ///
        double Coeff StdErr P_value CI_lower CI_upper ///
        using `box_params', replace

    foreach a of local age_levels {

        * Map age group numeric to English label directly via if
        local agelab ""
        if `a' == 1 local agelab "<=60 years"
        if `a' == 2 local agelab "61-70 years"
        if `a' == 3 local agelab "71-80 years"
        if `a' == 4 local agelab ">80 years"
        if `a' == 5 local agelab "All participants"

        foreach t of local type_levels {

            * Only compute for combinations with N >= 2
            count if age_group_ext == `a' & improve_type_en == "`t'"
            if r(N) >= 2 {
                quietly regress BAA_reduction ///
                    if age_group_ext == `a' & improve_type_en == "`t'"

                scalar b  = _b[_cons]
                scalar se = _se[_cons]
                scalar p  = 2 * ttail(e(df_r), abs(b/se))
                scalar tcrit = invttail(e(df_r), 0.025)
                scalar lb = b - tcrit * se
                scalar ub = b + tcrit * se

                local implab "`t'"
                post `memhold' ("`agelab'") ("`implab'") (b) (se) (p) (lb) (ub)
            }
        }
    }

    postclose `memhold'

    use `box_params', clear
    order AgeGroup ImprovementType Coeff StdErr P_value CI_lower CI_upper

    export delimited using ///
        "output/figures/counterfactual/baa_improvement_boxplot_params_by_age_and_type.csv", ///
        replace
restore
******************************************************


*------------------------------------------------------------
* Follow-up: age group descriptive statistics & other figures
*------------------------------------------------------------
display "=== Descriptive statistics of BAA improvement by age group ==="
tabstat BAA_reduction, by(age_group) stat(mean sd median min max n) format(%9.3f)

* ANOVA across age groups
oneway BAA_reduction age_group_order, tabulate
display "Age group F-test: F = " %6.2f r(F) ", P = " %6.4f r(p)

* Save age group analysis results
tempfile age_group_results
collapse (mean) baa_improvement_mean = BAA_reduction ///
         (sd)   baa_improvement_sd = BAA_reduction ///
         (count) n = BAA_reduction, by(age_group age_group_order)

sort age_group_order
gen se    = baa_improvement_sd / sqrt(n)
gen ci_lower = baa_improvement_mean - 1.96 * se
gen ci_upper = baa_improvement_mean + 1.96 * se

list age_group baa_improvement_mean se ci_lower ci_upper n, clean noobs
save "output/figures/counterfactual/age_group_baa_improvement_stats.dta", replace

* Reload valid data for subsequent figures
use "`valid_data'", clear
keep if !missing(BAA_reduction)

* Age group improvement bar chart
tempfile age_bar_data
collapse (mean) mean_BAA = BAA_reduction ///
         (sd)   sd_BAA   = BAA_reduction ///
         (count) n       = BAA_reduction, by(age_group age_group_order)

sort age_group_order
gen se = sd_BAA / sqrt(n)
gen lb = mean_BAA - 1.96 * se
gen ub = mean_BAA + 1.96 * se

graph bar mean_BAA, over(age_group_order, label(labsize(small))) ///
    bar(1, color("30 120 180")) ///
    blabel(bar, format(%4.3f)) ///
    ytitle("Mean BAA improvement (years)", size(small)) ///
    title("Mean BAA improvement by age group", size(medium)) ///
    note("Error bars indicate 95% confidence intervals", size(vsmall))
graph export "output/figures/counterfactual/mean_baa_improvement_bar_by_age.png", replace width(2000) height(1500)

* Reload valid data for subsequent figures
use "`valid_data'", clear
keep if !missing(BAA_reduction)

* 1. Improvement group bar chart
graph hbar (mean) BAA_reduction, over(improvement_type, sort(1)) ///
    title("Mean BAA reduction by improvement type", size(medium)) ///
    ytitle("BAA reduction (years)", size(small)) ///
    blabel(bar, format(%4.3f)) ///
    bargap(20) ///
    intensity(75)
graph export "output/figures/counterfactual/mean_baa_reduction_by_improvement_type.png", replace width(2000) height(1500)

* 2. Histogram of BAA reduction distribution
histogram BAA_reduction, ///
    title("Distribution of BAA reduction", size(medium)) ///
    xtitle("BAA reduction (years)", size(small)) ///
    ytitle("Frequency", size(small)) ///
    normal ///
    normopts(lcolor(red) lwidth(thick)) ///
    freq
graph export "output/figures/counterfactual/baa_reduction_histogram.png", replace width(2000) height(1500)

* 3. Box plot by year
graph box BAA_reduction, over(year) ///
    title("BAA reduction distribution by year", size(medium)) ///
    ytitle("BAA reduction (years)", size(small)) ///
    marker(1, msymbol(oh)) ///
    box(1, color(blue))
graph export "output/figures/counterfactual/baa_reduction_boxplot_by_year.png", replace width(2000) height(1500)

* 4. Box plot by improvement type
graph box BAA_reduction, over(improvement_type) ///
    title("BAA reduction distribution by improvement type", size(medium)) ///
    ytitle("BAA reduction (years)", size(small)) ///
    marker(1, msymbol(oh)) ///
    box(1, color(green)) ///
    intensity(60)
graph export "output/figures/counterfactual/baa_reduction_boxplot_by_improvement_type.png", replace width(2000) height(1500)

* 5. Scatter: EVI improvement vs BAA reduction
scatter BAA_reduction EVI_improvement, ///
    title("EVI improvement and BAA reduction", size(medium)) ///
    xtitle("EVI improvement", size(small)) ///
    ytitle("BAA reduction (years)", size(small)) ///
    mcolor(blue%50) ///
    || lfit BAA_reduction EVI_improvement, ///
    lcolor(red) lwidth(thick)
graph export "output/figures/counterfactual/evi_improvement_baa_reduction_scatter.png", replace width(2000) height(1500)

* 6. Scatter: O3 improvement vs BAA reduction
scatter BAA_reduction O3_improvement, ///
    title("O3 improvement and BAA reduction", size(medium)) ///
    xtitle("O3 improvement (ug/m3)", size(small)) ///
    ytitle("BAA reduction (years)", size(small)) ///
    mcolor(green%50) ///
    || lfit BAA_reduction O3_improvement, ///
    lcolor(red) lwidth(thick)
graph export "output/figures/counterfactual/o3_improvement_baa_reduction_scatter.png", replace width(2000) height(1500)

* 7. Pie chart of improvement type distribution
tempfile pie_data
use "`valid_data'", clear
keep if !missing(BAA_reduction)
contract improvement_type
gen proportion = _freq / _N * 100
graph pie proportion, over(improvement_type) ///
    title("Improvement type distribution", size(medium)) ///
    plabel(_all name, format(%2.1f) color(white)) ///
    legend(position(3))
graph export "output/figures/counterfactual/improvement_type_pie_chart.png", replace width(2000) height(1500)

* 8. Trend by year and improvement type
use "`valid_data'", clear
keep if !missing(BAA_reduction)
collapse (mean) BAA_reduction, by(year improvement_type)

twoway (connected BAA_reduction year if improvement_type == "both_improved", ///
           lcolor(red) mcolor(red)) ///
       (connected BAA_reduction year if improvement_type == "only_o3_improved", ///
           lcolor(blue) mcolor(blue)), ///
       title("BAA reduction trend by year and improvement type", size(medium)) ///
       xtitle("Year", size(small)) ///
       ytitle("Mean BAA reduction (years)", size(small)) ///
       legend(order(1 "Both improved" 2 "Only O3 improved") position(6) rows(1))
graph export "output/figures/counterfactual/baa_reduction_trend_by_year_and_type.png", replace width(2000) height(1500)

* Save final results dataset
use "`valid_data'", clear
keep if !missing(BAA_reduction)
keep id year age EVI_1KM_mean4y O3_1KM_mean4y counterfactual_EVI counterfactual_O3 ///
    actual_BAA_fixed counterfactual_BAA_fixed BAA_reduction baa_improved ///
    EVI_improvement O3_improvement improvement_type age_group age_group_order
save "output/tables/counterfactual_final_results.dta", replace

*------------------------------------------------------------
* Contribution of model coefficients to overall improvement
*------------------------------------------------------------
display "=== Contribution of model coefficients to improvement ==="

estimates replay double_exposure_model_main
matrix beta = e(b)
local evi_coef = beta[1,1]
local o3_coef  = beta[1,2]

display "EVI coefficient: `evi_coef'"
display "O3 coefficient: `o3_coef'"

* Calculate EVI contribution
summarize EVI_improvement
local evi_mean = r(mean)
display "Mean EVI improvement: `evi_mean'"
local evi_contribution = `evi_coef' * `evi_mean'
display "EVI contribution: `evi_coef' * `evi_mean' = `evi_contribution'"

* Calculate O3 contribution
summarize O3_improvement
local o3_mean = r(mean)
display "Mean O3 improvement: `o3_mean'"
local o3_contribution = `o3_coef' * `o3_mean'
display "O3 contribution: `o3_coef' * `o3_mean' = `o3_contribution'"

* Calculate total contribution
local total_contribution = `evi_contribution' + `o3_contribution'
display "Total contribution: `evi_contribution' + `o3_contribution' = `total_contribution'"

display ""
display "Average EVI improvement contribution: " %6.3f `evi_contribution' " years"
display "Average O3 improvement contribution: " %6.3f `o3_contribution' " years"
display "Total average contribution: " %6.3f `total_contribution' " years"

*------------------------------------------------------------
* Final summary
*------------------------------------------------------------
display ""
display "=== Analysis complete ==="
display "Model estimates saved as: double_exposure_model_main"
display "Counterfactual prediction results saved as: ../output/tables/counterfactual_final_results.dta"
display "Age group analysis results saved as: ../output/figures/counterfactual/age_group_baa_improvement_stats.dta"
display "Box plot parameters exported (English): ../output/figures/counterfactual/baa_improvement_boxplot_params_by_age_and_type.csv"
display ""
display "Exported figures saved to: ../output/figures/counterfactual/ folder"
display ""
display "Figure file list:"
display "1. baa_improvement_boxplot_by_age.png"
display "2. baa_improvement_boxplot_by_age_and_type.png (all English, includes all participants)"
display "3. mean_baa_improvement_bar_by_age.png"
display "4. mean_baa_reduction_by_improvement_type.png"
display "5. baa_reduction_histogram.png"
display "6. baa_reduction_boxplot_by_year.png"
display "7. baa_reduction_boxplot_by_improvement_type.png"
display "8. evi_improvement_baa_reduction_scatter.png"
display "9. o3_improvement_baa_reduction_scatter.png"
display "10. improvement_type_pie_chart.png"
display "11. baa_reduction_trend_by_year_and_type.png"
display ""
display "Analysis based on two key improvement groups: Both improved and Only O3 improved"
display "Age group analysis complete, includes 4 age groups: <=60 years, 61-70 years, 71-80 years, >80 years"
display "All English-required figures and parameter files exported as specified"

* Clean up temporary files (if still exist)
capture rm "`valid_data'"
capture rm "`age_group_results'"
capture rm "`age_bar_data'"
capture rm "`pie_data'"
