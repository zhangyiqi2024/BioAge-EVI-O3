set seed 20240101

use "DATA/analysis_ready/analysis_dataset.dta", clear

* Retain sample age 45-90
keep if age > 45 & age < 90

* Generate log-transformed household consumption
gen ln_hh_consumption_total = log(hh_consumption_total)

* Keep baseline data only (year=2011)
keep if year == 2011

* Create group labels
label define EVI_label 0 "Low EVI" 1 "High EVI"
label define O3_label 0 "Low O3" 1 "High O3"
label values traj15_EVI_1KM EVI_label
label values traj15_O3_1KM O3_label

* Check variable existence
describe BAA_PhenoAge BA_PhenoAge

* Climate variables principal component analysis
di "=== Climate variables principal component analysis ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* Save temporary data
save "temp_data.dta", replace

* Define outcome list and short names
local outcomes "BAA_PhenoAge BA_PhenoAge"
local shortnames "BAA BA"

* Loop analysis for each outcome
forvalues i = 1/2 {
    local outcome : word `i' of `outcomes'
    local outcome_name : word `i' of `shortnames'

    display "Analyzing outcome: `outcome' (`outcome_name')"

    * Load temporary data
    use "temp_data.dta", clear

    * Analysis 1: Marginal means for overall sample
    regress `outcome' i.traj15_EVI_1KM ///
        c.age i.sex i.education ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance ///
        i.smoking_current i.alcohol_current ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
        c.bmi i.hypertension i.dyslipidemia

    margins i.traj15_EVI_1KM, post
    parmest, saving("overall_margins_`outcome_name'.dta", replace)

    * Analysis 2: Marginal means for Low O3 group
    use "temp_data.dta", clear
    keep if traj15_O3_1KM == 0
    regress `outcome' i.traj15_EVI_1KM ///
        c.age i.sex i.education ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance ///
        i.smoking_current i.alcohol_current ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
        c.bmi i.hypertension i.dyslipidemia

    margins i.traj15_EVI_1KM, post
    parmest, saving("low_O3_margins_`outcome_name'.dta", replace)

    * Analysis 3: Marginal means for High O3 group
    use "temp_data.dta", clear
    keep if traj15_O3_1KM == 1
    regress `outcome' i.traj15_EVI_1KM ///
        c.age i.sex i.education ///
        c.ln_hh_consumption_total i.retirement i.pension_insurance ///
        i.smoking_current i.alcohol_current ///
        climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 ///
        c.bmi i.hypertension i.dyslipidemia

    margins i.traj15_EVI_1KM, post
    parmest, saving("high_O3_margins_`outcome_name'.dta", replace)

    * Combine all marginal means
    clear
    use "overall_margins_`outcome_name'.dta", clear
    gen group = "Overall"
    append using "low_O3_margins_`outcome_name'.dta"
    replace group = "Low O3 group" if mi(group)
    append using "high_O3_margins_`outcome_name'.dta"
    replace group = "High O3 group" if mi(group)

    * Extract EVI group information
    gen EVI_level = ""
    replace EVI_level = "Low EVI" if strpos(parm, "0b.traj15_EVI_1KM") == 0 & strpos(parm, "0.traj15_EVI_1KM") > 0
    replace EVI_level = "High EVI" if strpos(parm, "1.traj15_EVI_1KM") > 0
    keep if !mi(EVI_level)

    * Organize data
    keep group EVI_level estimate stderr min95 max95 p
    order group EVI_level estimate stderr min95 max95 p

    * Rename variables
    rename estimate marginal_mean
    rename stderr std_error
    rename min95 ci_lower
    rename max95 ci_upper
    rename p p_value

    * Add outcome identifier
    gen outcome_var = "`outcome_name'"

    * Create display group variable for plotting - in required order: Overall, Low O3, High O3
    gen display_group = ""
    replace display_group = "Overall" if group == "Overall"
    replace display_group = "Low O3" if group == "Low O3 group"
    replace display_group = "High O3" if group == "High O3 group"

    * Create group number - in required order: 1=Overall, 2=Low O3, 3=High O3
    gen group_num = .
    replace group_num = 1 if group == "Overall"
    replace group_num = 2 if group == "Low O3 group"
    replace group_num = 3 if group == "High O3 group"

    * Save results
    save "marginal_means_`outcome_name'.dta", replace

    display "Completed analysis for outcome: `outcome'"
}

* Combine results from both outcomes
clear
use "marginal_means_BAA.dta", clear
append using "marginal_means_BA.dta"

* Regenerate all figures to ensure files exist
display "=== Regenerating all figures to ensure files exist ==="

* Generate BAA horizontal box plot
preserve
keep if outcome_var == "BAA"
expand 100
bysort group EVI_level: gen simulated_value = marginal_mean + rnormal(0, std_error)

gen sort_order = .
replace sort_order = 1 if display_group == "Overall"
replace sort_order = 2 if display_group == "Low O3"
replace sort_order = 3 if display_group == "High O3"
label define group_label_order 1 "Overall" 2 "Low O3" 3 "High O3"
label values sort_order group_label_order

graph hbox simulated_value, ///
    over(EVI_level, gap(0)) ///
    over(sort_order, gap(50) label(labsize(medium))) ///
    ytitle("Adjusted BAA_PhenoAge", size(medium)) ///
    title("BAA_PhenoAge (horizontal box plot)", size(medium)) ///
    box(1, color(blue*0.7)) ///
    box(2, color(red*0.7)) ///
    intensity(70) ///
    legend(off)
graph save "output/figures/BAA_horizontal_boxplot.gph", replace
graph export "output/figures/BAA_horizontal_boxplot.png", as(png) width(3000) height(1200) replace
restore

* Generate BA horizontal box plot
preserve
keep if outcome_var == "BA"
expand 100
bysort group EVI_level: gen simulated_value = marginal_mean + rnormal(0, std_error)

gen sort_order = .
replace sort_order = 1 if display_group == "Overall"
replace sort_order = 2 if display_group == "Low O3"
replace sort_order = 3 if display_group == "High O3"
label define group_label_order 1 "Overall" 2 "Low O3" 3 "High O3"
label values sort_order group_label_order

graph hbox simulated_value, ///
    over(EVI_level, gap(0)) ///
    over(sort_order, gap(50) label(labsize(medium))) ///
    ytitle("Adjusted BA_PhenoAge", size(medium)) ///
    title("BA_PhenoAge (horizontal box plot)", size(medium)) ///
    box(1, color(blue*0.7)) ///
    box(2, color(red*0.7)) ///
    intensity(70) ///
    legend(off)
graph save "output/figures/BA_horizontal_boxplot.gph", replace
graph export "output/figures/BA_horizontal_boxplot.png", as(png) width(3000) height(1200) replace
restore

* Generate BAA vertical box plot
preserve
keep if outcome_var == "BAA"
expand 100
bysort group EVI_level: gen simulated_value = marginal_mean + rnormal(0, std_error)

gen sort_order = .
replace sort_order = 1 if display_group == "Overall"
replace sort_order = 2 if display_group == "Low O3"
replace sort_order = 3 if display_group == "High O3"
label define group_label_order 1 "Overall" 2 "Low O3" 3 "High O3"
label values sort_order group_label_order

graph box simulated_value, ///
    over(EVI_level, gap(0)) ///
    over(sort_order, gap(50) label(labsize(medium))) ///
    ytitle("Adjusted BAA_PhenoAge", size(medium)) ///
    title("BAA_PhenoAge (vertical box plot)", size(medium)) ///
    box(1, color(blue*0.7)) ///
    box(2, color(red*0.7)) ///
    intensity(70) ///
    legend(off)
graph save "output/figures/BAA_vertical_boxplot.gph", replace
graph export "output/figures/BAA_vertical_boxplot.png", as(png) width(3000) height(1200) replace
restore

* Generate BA vertical box plot
preserve
keep if outcome_var == "BA"
expand 100
bysort group EVI_level: gen simulated_value = marginal_mean + rnormal(0, std_error)

gen sort_order = .
replace sort_order = 1 if display_group == "Overall"
replace sort_order = 2 if display_group == "Low O3"
replace sort_order = 3 if display_group == "High O3"
label define group_label_order 1 "Overall" 2 "Low O3" 3 "High O3"
label values sort_order group_label_order

graph box simulated_value, ///
    over(EVI_level, gap(0)) ///
    over(sort_order, gap(50) label(labsize(medium))) ///
    ytitle("Adjusted BA_PhenoAge", size(medium)) ///
    title("BA_PhenoAge (vertical box plot)", size(medium)) ///
    box(1, color(blue*0.7)) ///
    box(2, color(red*0.7)) ///
    intensity(70) ///
    legend(off)
graph save "output/figures/BA_vertical_boxplot.gph", replace
graph export "output/figures/BA_vertical_boxplot.png", as(png) width(3000) height(1200) replace
restore

* Create combined horizontal box plot
graph combine "output/figures/BAA_horizontal_boxplot.gph" "output/figures/BA_horizontal_boxplot.gph", ///
    rows(1) ///
    title("Comparison of BAA_PhenoAge and BA_PhenoAge (horizontal box plots)", size(medium)) ///
    name("combined_horizontal_boxplot", replace)
graph save "output/figures/combined_horizontal_boxplot.gph", replace
graph export "output/figures/combined_horizontal_boxplot.png", as(png) width(3000) height(1200) replace

* Create combined vertical box plot
graph combine "output/figures/BAA_vertical_boxplot.gph" "output/figures/BA_vertical_boxplot.gph", ///
    rows(1) ///
    title("Comparison of BAA_PhenoAge and BA_PhenoAge (vertical box plots)", size(medium)) ///
    name("combined_vertical_boxplot", replace)
graph save "output/figures/combined_vertical_boxplot.gph", replace
graph export "output/figures/combined_vertical_boxplot.png", as(png) width(3000) height(1200) replace

* Create comprehensive comparison figure
graph combine "output/figures/combined_horizontal_boxplot.gph" "output/figures/combined_vertical_boxplot.gph", ///
    rows(2) ///
    title("Comprehensive comparison of BAA_PhenoAge and BA_PhenoAge", size(medium)) ///
    name("comprehensive_comparison", replace)
graph export "output/figures/comprehensive_comparison.png", as(png) width(3000) height(2000) replace

* Create detailed results table
preserve
collapse (mean) marginal_mean = marginal_mean (sd) std_error = std_error ///
    (mean) ci_lower = ci_lower (mean) ci_upper = ci_upper (mean) p_value = p_value, ///
    by(outcome_var group EVI_level)

* Add significance markers
gen significance = ""
replace significance = "***" if p_value < 0.001
replace significance = "**" if p_value >= 0.001 & p_value < 0.01
replace significance = "*" if p_value >= 0.01 & p_value < 0.05
replace significance = "" if p_value >= 0.05

* Format output
gen marginal_mean_formatted = string(marginal_mean, "%9.2f")
gen confidence_interval = "[" + string(ci_lower, "%9.2f") + ", " + string(ci_upper, "%9.2f") + "]"

keep outcome_var group EVI_level marginal_mean_formatted confidence_interval std_error significance
order outcome_var group EVI_level marginal_mean_formatted confidence_interval std_error significance

* Export to output directory
export excel using "output/tables/detailed_results_summary.xlsx", firstrow(variables) replace

* Also export as CSV for backup
export delimited using "output/tables/detailed_results_summary.csv", replace
restore

* Save final dataset
save "DATA/analysis_ready/final_analysis_dataset.dta", replace

* Output results summary
display "=== Comprehensive results: Adjusted baseline phenotypic age by EVI categories, overall and by O3 strata ==="
display ""
display "Multiple linear regression models were used to estimate marginal means of phenotypic age by EVI categories, adjusted for covariates."
display ""
list outcome_var group EVI_level marginal_mean ci_lower ci_upper p_value, sepby(outcome_var group) noobs

* Clean up temporary files
capture erase "temp_data.dta"
capture erase "overall_margins_BAA.dta"
capture erase "low_O3_margins_BAA.dta"
capture erase "high_O3_margins_BAA.dta"
capture erase "overall_margins_BA.dta"
capture erase "low_O3_margins_BA.dta"
capture erase "high_O3_margins_BA.dta"
capture erase "marginal_means_BAA.dta"
capture erase "marginal_means_BA.dta"

* Display completion info and file list
display ""
display "=== Comprehensive analysis completed ==="
display "All output files saved to ../output/ directory"
display ""
display "Output file list:"
display "1. ../output/figures/BAA_horizontal_boxplot.png - BAA_PhenoAge horizontal box plot"
display "2. ../output/figures/BA_horizontal_boxplot.png - BA_PhenoAge horizontal box plot"
display "3. ../output/figures/combined_horizontal_boxplot.png - Combined horizontal box plots"
display "4. ../output/figures/BAA_vertical_boxplot.png - BAA_PhenoAge vertical box plot"
display "5. ../output/figures/BA_vertical_boxplot.png - BA_PhenoAge vertical box plot"
display "6. ../output/figures/combined_vertical_boxplot.png - Combined vertical box plots"
display "7. ../output/figures/comprehensive_comparison.png - Comprehensive comparison"
display "8. ../output/tables/detailed_results_summary.xlsx - Detailed numerical results (Excel)"
display "9. ../output/tables/detailed_results_summary.csv - Detailed numerical results (CSV)"
display "10. ../DATA/analysis_ready/final_analysis_dataset.dta - Final analysis dataset"
