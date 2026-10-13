*******************************************************************************
// Analyze BAA_PhenoAge comparison across joint exposure groups using all data (2011 and 2015 waves)
* -----------------------------------------------------------------------------
* Step 4: Generate joint exposure group comparison plots - based on model-adjusted effect sizes
* -----------------------------------------------------------------------------

set seed 20240101

di as txt "=== Generating Joint Exposure Group Comparison Plots (Based on Model-Adjusted Effect Sizes) ==="

* Reload complete dataset
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
* Generate log-transformed household consumption variable
gen ln_hh_consumption_total = log(hh_consumption_total)
* Set graph style
set scheme s1color
graph set window fontface "Times New Roman"

* Principal Component Analysis of Climate Variables
di "=== Principal Component Analysis of Climate Variables ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* Method 1: Calculate adjusted marginal means using mixed-effects model
di as txt "Calculating model-adjusted marginal means..."

* 1KM scale - calculate adjusted means using mixed-effects model
mixed BAA_PhenoAge i.joint_exp_1KM                   ///
     c.age i.sex i.education                    ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
     i.smoking_current i.alcohol_current                           ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     c.bmi i.hypertension i.dyslipidemia               ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_1KM) vce(robust)

* Calculate marginal means and save results
margins joint_exp_1KM, saving(margins_1km, replace)
matrix margins_1km = r(table)

* 500m scale - calculate adjusted means using mixed-effects model
mixed BAA_PhenoAge i.joint_exp_500m                  ///
     c.age i.sex i.education                    ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance            ///
     i.smoking_current i.alcohol_current                           ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     c.bmi i.hypertension i.dyslipidemia               ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_500m) vce(robust)

* Calculate marginal means and save results
margins joint_exp_500m, saving(margins_500m, replace)
matrix margins_500m = r(table)

* Build plotting dataset
clear
set obs 8

* Create variables
generate scale = ""
generate group = ""
generate adj_mean = .
generate ci_low = .
generate ci_high = .
generate se = .
generate group_num = .

* Fill 1KM scale data
replace scale = "1KM" if _n <= 4
replace group_num = _n if _n <= 4

replace group = "Low EVI Low O3" if _n == 1
replace group = "Low EVI High O3" if _n == 2
replace group = "High EVI Low O3" if _n == 3
replace group = "High EVI High O3" if _n == 4

* Extract 1KM scale data from margins results
replace adj_mean = margins_1km[1,1] if _n == 1
replace se = margins_1km[2,1] if _n == 1
replace ci_low = margins_1km[5,1] if _n == 1
replace ci_high = margins_1km[6,1] if _n == 1

replace adj_mean = margins_1km[1,2] if _n == 2
replace se = margins_1km[2,2] if _n == 2
replace ci_low = margins_1km[5,2] if _n == 2
replace ci_high = margins_1km[6,2] if _n == 2

replace adj_mean = margins_1km[1,3] if _n == 3
replace se = margins_1km[2,3] if _n == 3
replace ci_low = margins_1km[5,3] if _n == 3
replace ci_high = margins_1km[6,3] if _n == 3

replace adj_mean = margins_1km[1,4] if _n == 4
replace se = margins_1km[2,4] if _n == 4
replace ci_low = margins_1km[5,4] if _n == 4
replace ci_high = margins_1km[6,4] if _n == 4

* Fill 500m scale data
replace scale = "500m" if _n > 4
replace group_num = _n - 4 if _n > 4

replace group = "Low EVI Low O3" if _n == 5
replace group = "Low EVI High O3" if _n == 6
replace group = "High EVI Low O3" if _n == 7
replace group = "High EVI High O3" if _n == 8

* Extract 500m scale data from margins results
replace adj_mean = margins_500m[1,1] if _n == 5
replace se = margins_500m[2,1] if _n == 5
replace ci_low = margins_500m[5,1] if _n == 5
replace ci_high = margins_500m[6,1] if _n == 5

replace adj_mean = margins_500m[1,2] if _n == 6
replace se = margins_500m[2,2] if _n == 6
replace ci_low = margins_500m[5,2] if _n == 6
replace ci_high = margins_500m[6,2] if _n == 6

replace adj_mean = margins_500m[1,3] if _n == 7
replace se = margins_500m[2,3] if _n == 7
replace ci_low = margins_500m[5,3] if _n == 7
replace ci_high = margins_500m[6,3] if _n == 7

replace adj_mean = margins_500m[1,4] if _n == 8
replace se = margins_500m[2,4] if _n == 8
replace ci_low = margins_500m[5,4] if _n == 8
replace ci_high = margins_500m[6,4] if _n == 8

* Generate plotting position variable
generate x_pos = group_num + (scale == "500m") * 0.4 - 0.2

* 1. Generate model-adjusted effect size comparison plot (two scales)
di as txt "Generating model-adjusted joint exposure group effect size comparison plot..."

twoway (bar adj_mean x_pos if scale == "1KM", barw(0.2) color(blue%70)) ///
       (bar adj_mean x_pos if scale == "500m", barw(0.2) color(red%70)) ///
       (rcap ci_high ci_low x_pos if scale == "1KM", lcolor(blue) lwidth(thin)) ///
       (rcap ci_high ci_low x_pos if scale == "500m", lcolor(red) lwidth(thin)), ///
       xlabel(1 "Low EVI Low O3" 2 "Low EVI High O3" 3 "High EVI Low O3" 4 "High EVI High O3", labsize(small)) ///
       ytitle("Adjusted BAA (years)", size(small) margin(medium)) ///
       xtitle("") ///
       legend(order(1 "1km scale" 2 "500m scale") rows(1) size(small)) ///
       xscale(range(0.5, 4.5)) ///
       yline(0, lcolor(black) lpattern(dash)) ///
       plotregion(margin(zero)) ///
       title("a", size(large) position(11) ring(0)) ///
       name(joint_exp, replace)

graph export "output/figures/joint_exposure_group_comparison_adjusted_means.png", replace width(2000) height(1200)
graph save joint_exp "output/figures/_joint_exp.gph", replace
di as txt "Joint exposure graph saved (PNG + GPH)"

* 2. Generate effect size forest plots for each scale separately
di as txt "Generating joint exposure group effect size forest plots..."

* 1KM scale forest plot
preserve
keep if scale == "1KM"

* Calculate effect values relative to reference group (Low EVI Low O3)
generate ref_mean = adj_mean[1]
generate effect_vs_ref = adj_mean - ref_mean
generate effect_se = sqrt(se^2 + se[1]^2)  // Simplified calculation

* Generate forest plot
twoway (scatter group_num effect_vs_ref, mcolor(blue) msymbol(O)) ///
       (rcap ci_high ci_low group_num, horizontal lcolor(blue)), ///
       ylabel(1 "Low EVI Low O3" 2 "Low EVI High O3" 3 "High EVI Low O3" 4 "High EVI High O3", angle(0)) ///
       xline(0, lcolor(black) lpattern(dash)) ///
       xtitle("Difference in BAA_PhenoAge compared to Low EVI Low O3 group", size(medsmall) margin(medium)) ///
       ytitle("") ///
       title("Forest Plot of Joint Exposure Effects on BAA_PhenoAge (1KM scale)", size(medium)) ///
       legend(off)

graph export "output/figures/forest_plot_1km_scale.png", replace width(1800) height(1000)
restore

* 500m scale forest plot
preserve
keep if scale == "500m"

* Calculate effect values relative to reference group (Low EVI Low O3)
generate ref_mean = adj_mean[1]
generate effect_vs_ref = adj_mean - ref_mean
generate effect_se = sqrt(se^2 + se[1]^2)

* Generate forest plot
twoway (scatter group_num effect_vs_ref, mcolor(red) msymbol(O)) ///
       (rcap ci_high ci_low group_num, horizontal lcolor(red)), ///
       ylabel(1 "Low EVI Low O3" 2 "Low EVI High O3" 3 "High EVI Low O3" 4 "High EVI High O3", angle(0)) ///
       xline(0, lcolor(black) lpattern(dash)) ///
       xtitle("Difference in BAA_PhenoAge compared to Low EVI Low O3 group", size(medsmall) margin(medium)) ///
       ytitle("") ///
       title("Forest Plot of Joint Exposure Effects on BAA_PhenoAge (500m scale)", size(medium)) ///
       legend(off)

graph export "output/figures/forest_plot_500m_scale.png", replace width(1800) height(1000)
restore

* 3. Generate effect size table output
di as txt "Outputting joint exposure group effect size table..."

* Save effect size results
preserve
keep scale group adj_mean se ci_low ci_high
gen stars = ""
replace stars = "*" if (ci_low > 0 & ci_high > 0) | (ci_low < 0 & ci_high < 0)
replace stars = "**" if abs(adj_mean/se) > 2.58  // approximate p<0.01

format adj_mean se ci_low ci_high %8.3f

list scale group adj_mean se ci_low ci_high stars, sepby(scale) noobs

* Output to file
export delimited using "output/tables/joint_exposure_group_effects.csv", replace
restore

* 4. Generate between-group difference comparison plots
di as txt "Generating joint exposure group difference comparison plots..."

* Re-run model for formal between-group comparisons
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90

* Generate log-transformed household consumption variable
gen ln_hh_consumption_total = log(hh_consumption_total)

* Principal Component Analysis of Climate Variables
di "=== Principal Component Analysis of Climate Variables ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* 1KM scale between-group comparison
mixed BAA_PhenoAge i.joint_exp_1KM                   ///
     age i.sex i.education                    ///
     ln_hh_consumption_total i.retirement i.pension_insurance            ///
     i.smoking_current i.alcohol_current                           ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     bmi i.hypertension i.dyslipidemia               ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_1KM) vce(robust)

* Conduct all pairwise between-group comparisons
margins joint_exp_1KM, pwcompare(effects) mcompare(bonferroni)

* Save comparison results
matrix comp_1km = r(table)

* Output between-group comparison results
di as txt "=== 1KM Scale Group Comparison Results ==="
margins joint_exp_1KM, pwcompare

* 500m scale between-group comparison
mixed BAA_PhenoAge i.joint_exp_500m                  ///
     age i.sex i.education                    ///
     ln_hh_consumption_total i.retirement i.pension_insurance            ///
     i.smoking_current i.alcohol_current                           ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     bmi i.hypertension i.dyslipidemia               ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_500m) vce(robust)

* Conduct all pairwise between-group comparisons
margins joint_exp_500m, pwcompare(effects) mcompare(bonferroni)

* Output between-group comparison results
di as txt "=== 500m Scale Group Comparison Results ==="
margins joint_exp_500m, pwcompare

di as txt "=== Joint Exposure Group Comparison Analysis Completed! ==="
di as txt "Generated plots:"
di as txt "1. joint_exposure_group_comparison_adjusted_means.png - Comparison of adjusted means across two scales"
di as txt "2. forest_plot_1km_scale.png - Forest plot of effect sizes at 1KM scale"
di as txt "3. forest_plot_500m_scale.png - Forest plot of effect sizes at 500m scale"
di as txt "4. joint_exposure_group_effects.csv - Detailed effect size table"

* Clean up temporary files
capture erase "margins_1km.dta"
capture erase "margins_500m.dta"
