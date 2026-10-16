* -----------------------------------------------------------------------------
* Complete code: Run all regression models and generate exposure effect summary table
* (including dual-exposure models)
* -----------------------------------------------------------------------------
set seed 20240101

* Reload data (working directory must be the project root: the folder
* containing CODE/, DATA/ and output/ - see CODE/00_master_run.do)
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)

* Sort by individual and year
sort id_num year

* Climate variable principal component analysis
di "=== Climate Variable PCA ==="
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* -----------------------------------------------------------------------------
* Step 1: Run all regression models
* -----------------------------------------------------------------------------

di as txt "=== Starting all regression models ==="

*********************** Model 1-1: EVI_1KM BAA_PhenoAge (mixed effects)
*---- Mixed effects model (weighted)
di as txt "=== Model 1-1: EVI_1KM BAA_PhenoAge (weighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_1KM                  ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_EVI_1KM) vce(robust)
estimates store model1_1_pwe

*---- Mixed effects model (unweighted)
di as txt "=== Model 1-1: EVI_1KM BAA_PhenoAge (unweighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_1KM                  ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     vce(robust)
estimates store model1_1_0

*---- Marginal structural model (MSM)
di as txt "=== Model 1-1: EVI_1KM BAA_PhenoAge (marginal structural model) ==="
regress BAA_PhenoAge i.traj15_EVI_1KM                ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     [pweight=ipw_EVI_1KM],                          ///
     vce(cluster id_num)
estimates store model1_1_msm

*********************** Model 1-2: EVI_500m BAA_PhenoAge (mixed effects)
*---- Mixed effects model (weighted)
di as txt "=== Model 1-2: EVI_500m BAA_PhenoAge (weighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_500m                 ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_EVI_500m) vce(robust)
estimates store model1_2_pwe

*---- Mixed effects model (unweighted)
di as txt "=== Model 1-2: EVI_500m BAA_PhenoAge (unweighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_500m                 ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     vce(robust)
estimates store model1_2_0

*---- Marginal structural model (MSM)
di as txt "=== Model 1-2: EVI_500m BAA_PhenoAge (marginal structural model) ==="
regress BAA_PhenoAge i.traj15_EVI_500m               ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     [pweight=ipw_EVI_500m],                         ///
     vce(cluster id_num)
estimates store model1_2_msm

*********************** Model 2-1: O3_1KM BAA_PhenoAge (mixed effects)
*---- Mixed effects model (weighted)
di as txt "=== Model 2-1: O3_1KM BAA_PhenoAge (weighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_O3_1KM                   ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_O3_1KM) vce(robust)
estimates store model2_1_pwe

*---- Mixed effects model (unweighted)
di as txt "=== Model 2-1: O3_1KM BAA_PhenoAge (unweighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_O3_1KM                   ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     vce(robust)
estimates store model2_1_0

*---- Marginal structural model (MSM)
di as txt "=== Model 2-1: O3_1KM BAA_PhenoAge (marginal structural model) ==="
regress BAA_PhenoAge i.traj15_O3_1KM                 ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     [pweight=ipw_O3_1KM],                           ///
     vce(cluster id_num)
estimates store model2_1_msm

*********************** Model 2-2: O3_500m BAA_PhenoAge (mixed effects)
*---- Mixed effects model (weighted)
di as txt "=== Model 2-2: O3_500m BAA_PhenoAge (weighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_O3_500m                  ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_O3_500m) vce(robust)
estimates store model2_2_pwe

*---- Mixed effects model (unweighted)
di as txt "=== Model 2-2: O3_500m BAA_PhenoAge (unweighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_O3_500m                  ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     vce(robust)
estimates store model2_2_0

*---- Marginal structural model (MSM)
di as txt "=== Model 2-2: O3_500m BAA_PhenoAge (marginal structural model) ==="
regress BAA_PhenoAge i.traj15_O3_500m                ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     [pweight=ipw_O3_500m],                          ///
     vce(cluster id_num)
estimates store model2_2_msm

* -----------------------------------------------------------------------------
* Dual-exposure regression models
* -----------------------------------------------------------------------------

*********************** Model 3-1: EVI_1KM O3_1KM (dual-exposure)
*---- Mixed effects model (weighted)
di as txt "=== Model 3-1: EVI_1KM + O3_1KM BAA_PhenoAge (weighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_1KM i.traj15_O3_1KM  ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_1KM) vce(robust)
estimates store model3_1_pwe

*---- Mixed effects model (unweighted)
di as txt "=== Model 3-1: EVI_1KM + O3_1KM BAA_PhenoAge (unweighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_1KM i.traj15_O3_1KM  ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     vce(robust)
estimates store model3_1_0

*---- Marginal structural model (MSM)
di as txt "=== Model 3-1: EVI_1KM + O3_1KM BAA_PhenoAge (marginal structural model) ==="
regress BAA_PhenoAge i.traj15_EVI_1KM i.traj15_O3_1KM ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     [pweight=ipw_joint_1KM],                        ///
     vce(cluster id_num)
estimates store model3_1_msm

*********************** Model 3-2: EVI_500m O3_500m (dual-exposure model)
*---- Mixed effects model (weighted)
di as txt "=== Model 3-2: EVI_500m + O3_500m BAA_PhenoAge (weighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_500m i.traj15_O3_500m ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     pweight(ipw_joint_500m) vce(robust)
estimates store model3_2_pwe

*---- Mixed effects model (unweighted)
di as txt "=== Model 3-2: EVI_500m + O3_500m BAA_PhenoAge (unweighted mixed effects) ==="
mixed BAA_PhenoAge i.traj15_EVI_500m i.traj15_O3_500m ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
	 climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
	 c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     || id_num:,                                     ///
     vce(robust)
estimates store model3_2_0

*---- Marginal structural model (MSM)
di as txt "=== Model 3-2: EVI_500m + O3_500m BAA_PhenoAge (marginal structural model) ==="
regress BAA_PhenoAge i.traj15_EVI_500m i.traj15_O3_500m ///
     c.age i.sex i.education                        ///
     c.ln_hh_consumption_total i.retirement i.pension_insurance  ///
     i.smoking_current i.alcohol_current            ///
     climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5  ///
     c.bmi i.hypertension i.dyslipidemia            ///
     i.year                                          ///
     [pweight=ipw_joint_500m],                       ///
     vce(cluster id_num)
estimates store model3_2_msm

di as green "=== All regression models completed ==="

* -----------------------------------------------------------------------------
* Step 2: Use lincom to extract exposure variable coefficients for all models
* -----------------------------------------------------------------------------

di as txt "=== Extracting exposure variable coefficients for all models ==="

* Create empty dataset to store results - now 24 models (12 single-exposure + 12 dual-exposure)
clear
set obs 24

* Create variables
gen model_id = _n
gen model_description = ""
gen exposure_variable = ""
gen coefficient = .
gen std_error = .
gen p_value = .
gen sample_size = .
gen significance = ""

* Extract coefficients in specified order

* ==================== Single-exposure models ====================

* EVI 1KM WME
estimates restore model1_1_pwe
lincom 1.traj15_EVI_1KM
replace coefficient = r(estimate) in 1
replace std_error = r(se) in 1
replace p_value = r(p) in 1
replace sample_size = e(N) in 1
replace model_description = "EVI 1KM WME" in 1
replace exposure_variable = "traj15_EVI_1KM" in 1

* EVI 1KM UME
estimates restore model1_1_0
lincom 1.traj15_EVI_1KM
replace coefficient = r(estimate) in 2
replace std_error = r(se) in 2
replace p_value = r(p) in 2
replace sample_size = e(N) in 2
replace model_description = "EVI 1KM UME" in 2
replace exposure_variable = "traj15_EVI_1KM" in 2

* EVI 1KM MSM
estimates restore model1_1_msm
lincom 1.traj15_EVI_1KM
replace coefficient = r(estimate) in 3
replace std_error = r(se) in 3
replace p_value = r(p) in 3
replace sample_size = e(N) in 3
replace model_description = "EVI 1KM MSM" in 3
replace exposure_variable = "traj15_EVI_1KM" in 3

* EVI 500M WME
estimates restore model1_2_pwe
lincom 1.traj15_EVI_500m
replace coefficient = r(estimate) in 4
replace std_error = r(se) in 4
replace p_value = r(p) in 4
replace sample_size = e(N) in 4
replace model_description = "EVI 500M WME" in 4
replace exposure_variable = "traj15_EVI_500m" in 4

* EVI 500M UME
estimates restore model1_2_0
lincom 1.traj15_EVI_500m
replace coefficient = r(estimate) in 5
replace std_error = r(se) in 5
replace p_value = r(p) in 5
replace sample_size = e(N) in 5
replace model_description = "EVI 500M UME" in 5
replace exposure_variable = "traj15_EVI_500m" in 5

* EVI 500M MSM
estimates restore model1_2_msm
lincom 1.traj15_EVI_500m
replace coefficient = r(estimate) in 6
replace std_error = r(se) in 6
replace p_value = r(p) in 6
replace sample_size = e(N) in 6
replace model_description = "EVI 500M MSM" in 6
replace exposure_variable = "traj15_EVI_500m" in 6

* O3 1KM WME
estimates restore model2_1_pwe
lincom 1.traj15_O3_1KM
replace coefficient = r(estimate) in 7
replace std_error = r(se) in 7
replace p_value = r(p) in 7
replace sample_size = e(N) in 7
replace model_description = "O3 1KM WME" in 7
replace exposure_variable = "traj15_O3_1KM" in 7

* O3 1KM UME
estimates restore model2_1_0
lincom 1.traj15_O3_1KM
replace coefficient = r(estimate) in 8
replace std_error = r(se) in 8
replace p_value = r(p) in 8
replace sample_size = e(N) in 8
replace model_description = "O3 1KM UME" in 8
replace exposure_variable = "traj15_O3_1KM" in 8

* O3 1KM MSM
estimates restore model2_1_msm
lincom 1.traj15_O3_1KM
replace coefficient = r(estimate) in 9
replace std_error = r(se) in 9
replace p_value = r(p) in 9
replace sample_size = e(N) in 9
replace model_description = "O3 1KM MSM" in 9
replace exposure_variable = "traj15_O3_1KM" in 9

* O3 500M WME
estimates restore model2_2_pwe
lincom 1.traj15_O3_500m
replace coefficient = r(estimate) in 10
replace std_error = r(se) in 10
replace p_value = r(p) in 10
replace sample_size = e(N) in 10
replace model_description = "O3 500M WME" in 10
replace exposure_variable = "traj15_O3_500m" in 10

* O3 500M UME
estimates restore model2_2_0
lincom 1.traj15_O3_500m
replace coefficient = r(estimate) in 11
replace std_error = r(se) in 11
replace p_value = r(p) in 11
replace sample_size = e(N) in 11
replace model_description = "O3 500M UME" in 11
replace exposure_variable = "traj15_O3_500m" in 11

* O3 500M MSM
estimates restore model2_2_msm
lincom 1.traj15_O3_500m
replace coefficient = r(estimate) in 12
replace std_error = r(se) in 12
replace p_value = r(p) in 12
replace sample_size = e(N) in 12
replace model_description = "O3 500M MSM" in 12
replace exposure_variable = "traj15_O3_500m" in 12

* ==================== Dual-exposure models ====================

* Joint EVI 1KM WME
estimates restore model3_1_pwe
lincom 1.traj15_EVI_1KM
replace coefficient = r(estimate) in 13
replace std_error = r(se) in 13
replace p_value = r(p) in 13
replace sample_size = e(N) in 13
replace model_description = "Joint EVI 1KM WME" in 13
replace exposure_variable = "traj15_EVI_1KM_joint" in 13

* Joint O3 1KM WME
estimates restore model3_1_pwe
lincom 1.traj15_O3_1KM
replace coefficient = r(estimate) in 14
replace std_error = r(se) in 14
replace p_value = r(p) in 14
replace sample_size = e(N) in 14
replace model_description = "Joint O3 1KM WME" in 14
replace exposure_variable = "traj15_O3_1KM_joint" in 14

* Joint EVI 1KM UME
estimates restore model3_1_0
lincom 1.traj15_EVI_1KM
replace coefficient = r(estimate) in 15
replace std_error = r(se) in 15
replace p_value = r(p) in 15
replace sample_size = e(N) in 15
replace model_description = "Joint EVI 1KM UME" in 15
replace exposure_variable = "traj15_EVI_1KM_joint" in 15

* Joint O3 1KM UME
estimates restore model3_1_0
lincom 1.traj15_O3_1KM
replace coefficient = r(estimate) in 16
replace std_error = r(se) in 16
replace p_value = r(p) in 16
replace sample_size = e(N) in 16
replace model_description = "Joint O3 1KM UME" in 16
replace exposure_variable = "traj15_O3_1KM_joint" in 16

* Joint EVI 1KM MSM
estimates restore model3_1_msm
lincom 1.traj15_EVI_1KM
replace coefficient = r(estimate) in 17
replace std_error = r(se) in 17
replace p_value = r(p) in 17
replace sample_size = e(N) in 17
replace model_description = "Joint EVI 1KM MSM" in 17
replace exposure_variable = "traj15_EVI_1KM_joint" in 17

* Joint O3 1KM MSM
estimates restore model3_1_msm
lincom 1.traj15_O3_1KM
replace coefficient = r(estimate) in 18
replace std_error = r(se) in 18
replace p_value = r(p) in 18
replace sample_size = e(N) in 18
replace model_description = "Joint O3 1KM MSM" in 18
replace exposure_variable = "traj15_O3_1KM_joint" in 18

* Joint EVI 500M WME
estimates restore model3_2_pwe
lincom 1.traj15_EVI_500m
replace coefficient = r(estimate) in 19
replace std_error = r(se) in 19
replace p_value = r(p) in 19
replace sample_size = e(N) in 19
replace model_description = "Joint EVI 500M WME" in 19
replace exposure_variable = "traj15_EVI_500m_joint" in 19

* Joint O3 500M WME
estimates restore model3_2_pwe
lincom 1.traj15_O3_500m
replace coefficient = r(estimate) in 20
replace std_error = r(se) in 20
replace p_value = r(p) in 20
replace sample_size = e(N) in 20
replace model_description = "Joint O3 500M WME" in 20
replace exposure_variable = "traj15_O3_500m_joint" in 20

* Joint EVI 500M UME
estimates restore model3_2_0
lincom 1.traj15_EVI_500m
replace coefficient = r(estimate) in 21
replace std_error = r(se) in 21
replace p_value = r(p) in 21
replace sample_size = e(N) in 21
replace model_description = "Joint EVI 500M UME" in 21
replace exposure_variable = "traj15_EVI_500m_joint" in 21

* Joint O3 500M UME
estimates restore model3_2_0
lincom 1.traj15_O3_500m
replace coefficient = r(estimate) in 22
replace std_error = r(se) in 22
replace p_value = r(p) in 22
replace sample_size = e(N) in 22
replace model_description = "Joint O3 500M UME" in 22
replace exposure_variable = "traj15_O3_500m_joint" in 22

* Joint EVI 500M MSM
estimates restore model3_2_msm
lincom 1.traj15_EVI_500m
replace coefficient = r(estimate) in 23
replace std_error = r(se) in 23
replace p_value = r(p) in 23
replace sample_size = e(N) in 23
replace model_description = "Joint EVI 500M MSM" in 23
replace exposure_variable = "traj15_EVI_500m_joint" in 23

* Joint O3 500M MSM
estimates restore model3_2_msm
lincom 1.traj15_O3_500m
replace coefficient = r(estimate) in 24
replace std_error = r(se) in 24
replace p_value = r(p) in 24
replace sample_size = e(N) in 24
replace model_description = "Joint O3 500M MSM" in 24
replace exposure_variable = "traj15_O3_500m_joint" in 24

di as green "=== All model coefficients extracted successfully ==="

* -----------------------------------------------------------------------------
* Step 3: Generate final results table and forest plot
* -----------------------------------------------------------------------------

* Calculate confidence intervals and significance
gen CI_lower = coefficient - 1.96 * std_error
gen CI_upper = coefficient + 1.96 * std_error

replace significance = ""
replace significance = "*" if p_value < 0.05
replace significance = "**" if p_value < 0.01
replace significance = "***" if p_value < 0.001

* Add classification variables
gen exposure_type = ""
replace exposure_type = "Single EVI" if inlist(model_id, 1, 2, 3, 4, 5, 6)
replace exposure_type = "Single O3" if inlist(model_id, 7, 8, 9, 10, 11, 12)
replace exposure_type = "Joint EVI" if strpos(model_description, "Joint EVI")
replace exposure_type = "Joint O3" if strpos(model_description, "Joint O3")

gen buffer_distance = ""
replace buffer_distance = "1KM" if inlist(model_id, 1, 2, 3, 7, 8, 9, 13, 14, 15, 16, 17, 18)
replace buffer_distance = "500M" if inlist(model_id, 4, 5, 6, 10, 11, 12, 19, 20, 21, 22, 23, 24)

gen analysis_method = ""
replace analysis_method = "WME" if strpos(model_description, "WME")
replace analysis_method = "UME" if strpos(model_description, "UME")
replace analysis_method = "MSM" if strpos(model_description, "MSM")

gen outcome_type = "BAA"

* Format output
gen coefficient_CI = string(coefficient, "%9.3f") + significance + " (" + string(CI_lower, "%9.3f") + ", " + string(CI_upper, "%9.3f") + ")"
gen SE = string(std_error, "%9.3f")

* Reorder variables
order model_id model_description exposure_type buffer_distance analysis_method coefficient_CI SE coefficient std_error p_value

* Output table
di as green "=== Summary Table of Exposure Effects (Ordered) ==="
list model_id model_description coefficient_CI SE if model_id <= 24, noobs clean

* Export to Excel
preserve
    keep model_description exposure_type buffer_distance analysis_method coefficient_CI SE coefficient std_error p_value sample_size
    export excel using "output/tables/Ordered_Exposure_Effects_Summary_Table.xlsx", firstrow(variables) replace
restore

* -----------------------------------------------------------------------------
* Step 4: Create forest plot (ordered)
* -----------------------------------------------------------------------------

* Create y-axis positions for forest plot
gen y_position = _n

* Create model labels for forest plot (ordered)
gen model_label = model_description

* Set graph scheme
set scheme s1color

* Create forest plot - using exact specified order
twoway (rcap CI_lower CI_upper y_position, horizontal lcolor(navy) lwidth(0.8)) ///
       (scatter y_position coefficient, mcolor(maroon) msymbol(square) msize(medium)) ///
       , xline(0, lpattern(dash) lcolor(gray)) ///
         title("Associations Between Environmental Exposures and Biological Age Acceleration", size(medium)) ///
         xtitle("Coefficient (95% CI)", size(small)) ///
         ytitle("") ///
         xlabel(-3 -2 -1 0 1 2 3, labsize(small)) ///
         ylabel(1 "EVI 1KM WME"    2 "EVI 1KM UME"     3 "EVI 1KM MSM" ///
                4 "EVI 500M WME"   5 "EVI 500M UME"    6 "EVI 500M MSM" ///
                7 "O3 1KM WME"     8 "O3 1KM UME"      9 "O3 1KM MSM" ///
                10 "O3 500M WME"   11 "O3 500M UME"    12 "O3 500M MSM" ///
                13 "Joint EVI 1KM WME"   14 "Joint O3 1KM WME" ///
                15 "Joint EVI 1KM UME"   16 "Joint O3 1KM UME" ///
                17 "Joint EVI 1KM MSM"   18 "Joint O3 1KM MSM" ///
                19 "Joint EVI 500M WME"  20 "Joint O3 500M WME" ///
                21 "Joint EVI 500M UME"  22 "Joint O3 500M UME" ///
                23 "Joint EVI 500M MSM"  24 "Joint O3 500M MSM", ///
                angle(0) labsize(vsmall)) ///
         legend(off) ///
         graphregion(color(white)) ///
         plotregion(color(white)) ///
         yscale(reverse) ///
         name(ordered_forest_plot, replace)

* Save forest plot
graph export "output/figures/Ordered_Combined_Forest_Plot.png", width(2000) replace
graph export "output/figures/Ordered_Combined_Forest_Plot.pdf", replace

* Create detailed grouped summary table
preserve
    gen category = exposure_type + " " + buffer_distance + " " + analysis_method
    keep category model_description coefficient_CI SE coefficient std_error p_value sample_size
    order category model_description coefficient_CI SE
    export excel using "output/tables/Detailed_Ordered_Results_Table.xlsx", firstrow(variables) replace
restore

di as green "=== Analysis Complete ==="
di as green "=== Ordered forest plot saved as: ../output/figures/Ordered_Combined_Forest_Plot.png and .pdf ==="
di as green "=== Ordered summary table saved as: ../output/tables/Ordered_Exposure_Effects_Summary_Table.xlsx ==="
di as green "=== Detailed ordered table saved as: ../output/tables/Detailed_Ordered_Results_Table.xlsx ==="

* Create grouped statistics
di as green "=== Results by Exposure Type ==="
table exposure_type analysis_method, statistic(mean coefficient) statistic(mean sample_size)

di as green "=== Results by Buffer Distance ==="
table buffer_distance exposure_type, statistic(mean coefficient) statistic(mean sample_size)
