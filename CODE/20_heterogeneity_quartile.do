* =============================================================================
* Complete EVI and O3 Dual-Exposure Model Analysis - Quartile Score Version
* Includes: Data Preparation, Model Analysis, Interaction, Stratified Analysis, Heterogeneity Test, Result Export
* =============================================================================

set seed 20240101

di as txt "=== EVI and O3 Dual-Exposure Model Analysis - Complete Code ==="

* 0. Set Working Directory and Graph Style ---------------------------------------------------
* (working directory must be the project root: the folder containing CODE/,
*  DATA/ and output/ - see CODE/00_master_run.do)
set scheme s1color
graph set window fontface "Times New Roman"

* 1. Data Preparation and Preprocessing ---------------------------------------------------------
di as txt "=== 1. Data Preparation and Preprocessing ==="

use "DATA/analysis_ready/analysis_dataset.dta", clear

* Check key variables
di "=== 关键变量检查 ==="
sum EVI_1KM_mean4y O3_1KM_mean4y BAA_PhenoAge

* 限制age范围
keep if age > 45 & age < 90
di "age筛选后Sample size: " _N

* Log household consumption
gen ln_hh_consumption_total = log(hh_consumption_total)

* age_binary
gen age_binary = 0 if age <= 65
replace age_binary = 1 if age > 65 & age < .

* 2. Calculate quartile scores --------------------------------------------------------
di as txt "=== 2. Calculate quartile scores ==="

* Calculate EVI and O3 IQR
sum EVI_1KM_mean4y, detail
local evi_iqr = r(p75) - r(p25)
di "EVI IQR: `evi_iqr'"

sum O3_1KM_mean4y, detail  
local o3_iqr = r(p75) - r(p25)
di "O3 IQR: `o3_iqr'"

* Generate quartile variables
xtile EVI_quartile = EVI_1KM_mean4y, nq(4)
xtile O3_quartile = O3_1KM_mean4y, nq(4)

* Generate score variables (0-3 pts)
* EVI: higher is better, Q4 (highest) = 3pts, Q1 (lowest) = 0pts
gen EVI_score = EVI_quartile - 1

* O3: lower is better, Q1 (lowest) = 3pts, Q4 (highest) = 0pts  
gen O3_score = 3 - (O3_quartile - 1)

* Calculate total score (0-6 pts)
gen total_exposure_score = EVI_score + O3_score

* Check score distribution
tab EVI_score
tab O3_score  
tab total_exposure_score

* 3. Climate PCA -------------------------------------------------------
di as txt "=== 3. Climate PCA ==="

pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

* Define control variable list
local controls "c.age i.sex i.education c.ln_hh_consumption_total i.retirement i.pension_insurance i.smoking_current i.alcohol_current climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 c.bmi i.hypertension i.dyslipidemia i.year"

* 4. Full sample mixed exposure model ------------------------------------------------------
di as txt "=== 4. Full sample mixed exposure model ==="

mixed BAA_PhenoAge c.total_exposure_score `controls' || id_num:

* Extract coefficient and standard error
matrix b = e(b)
matrix V = e(V)

local score_coef_full = b[1,1]
local score_se_full = sqrt(V[1,1])
local score_z_full = `score_coef_full' / `score_se_full'
local score_p_full = 2 * (1 - normal(abs(`score_z_full')))
local n_full = e(N)

di "Full sample mixed exposure model结果:"
di "Total score coefficient: " %8.4f `score_coef_full' " (SE: " %8.4f `score_se_full' ", p: " %6.4f `score_p_full' ")"
di "Sample size: `n_full'"

* 5. Interaction Analysis --------------------------------------------------------------
di as txt "=== 5. Interaction Analysis ==="

* 5.1 sex × 混合暴露得分
di "--- sex交互作用 ---"
mixed BAA_PhenoAge c.total_exposure_score##i.sex `controls' || id_num:
estimates store model_gender_interaction
testparm i.sex#c.total_exposure_score
local p_interaction_gender = r(p)
di "sex交互作用p值: `p_interaction_gender'"

* 5.2 age分组 × 混合暴露得分
di "--- age交互作用 ---"
mixed BAA_PhenoAge c.total_exposure_score##i.age_binary `controls' || id_num:
estimates store model_age_interaction
testparm i.age_binary#c.total_exposure_score
local p_interaction_age = r(p)
di "age交互作用p值: `p_interaction_age'"

* 5.3 Region x Mixed Exposure Score
di "--- Region interaction ---"
mixed BAA_PhenoAge c.total_exposure_score##i.region_class2 `controls' || id_num:
estimates store model_region_interaction
testparm i.region_class2#c.total_exposure_score
local p_interaction_region = r(p)
di "Region interactionp值: `p_interaction_region'"

* 6. Stratified analysis ----------------------------------------------------------------
di as txt "=== 6. Stratified analysis ==="

* 6.1 sex分层
di "--- sex分层 ---"

preserve
keep if sex == 1
mixed BAA_PhenoAge c.total_exposure_score `controls' || id_num:
matrix b_male = e(b)
matrix V_male = e(V)
local score_coef_male = b_male[1,1]
local score_se_male = sqrt(V_male[1,1])
local n_male = e(N)
di "MaleSample size: `n_male', 系数: `score_coef_male', 标准误: `score_se_male'"
restore

preserve
keep if sex == 0
mixed BAA_PhenoAge c.total_exposure_score `controls' || id_num:
matrix b_female = e(b)
matrix V_female = e(V)
local score_coef_female = b_female[1,1]
local score_se_female = sqrt(V_female[1,1])
local n_female = e(N)
di "FemaleSample size: `n_female', 系数: `score_coef_female', 标准误: `score_se_female'"
restore

* 6.2 age分层
di "--- age分层 ---"

preserve
keep if age_binary == 0
mixed BAA_PhenoAge c.total_exposure_score i.sex `controls' || id_num:
matrix b_young = e(b)
matrix V_young = e(V)
local score_coef_young = b_young[1,1]
local score_se_young = sqrt(V_young[1,1])
local n_young = e(N)
di "<=65 yrsSample size: `n_young', 系数: `score_coef_young', 标准误: `score_se_young'"
restore

preserve
keep if age_binary == 1
mixed BAA_PhenoAge c.total_exposure_score i.sex `controls' || id_num:
matrix b_old = e(b)
matrix V_old = e(V)
local score_coef_old = b_old[1,1]
local score_se_old = sqrt(V_old[1,1])
local n_old = e(N)
di ">65 yrsSample size: `n_old', 系数: `score_coef_old', 标准误: `score_se_old'"
restore

* 6.3 Region stratification
di "--- Region stratification ---"

preserve
keep if region_class2 == 1
mixed BAA_PhenoAge c.total_exposure_score `controls' || id_num:
matrix b_north = e(b)
matrix V_north = e(V)
local score_coef_north = b_north[1,1]
local score_se_north = sqrt(V_north[1,1])
local n_north = e(N)
di "NorthSample size: `n_north', 系数: `score_coef_north', 标准误: `score_se_north'"
restore

preserve
keep if region_class2 == 2
mixed BAA_PhenoAge c.total_exposure_score `controls' || id_num:
matrix b_south = e(b)
matrix V_south = e(V)
local score_coef_south = b_south[1,1]
local score_se_south = sqrt(V_south[1,1])
local n_south = e(N)
di "SouthSample size: `n_south', 系数: `score_coef_south', 标准误: `score_se_south'"
restore

* 7. Formal test of subgroup differences ---------------------------------------------------
di as txt "=== 7. Formal test of subgroup differences ==="

* 7.1 Likelihood ratio test
di "--- Likelihood ratio test ---"

* Model without interactions (constrained)
mixed BAA_PhenoAge c.total_exposure_score i.sex i.age_binary i.region_class2 `controls' || id_num:
estimates store model_constrained

* Model with interactions (full)
mixed BAA_PhenoAge c.total_exposure_score##i.sex c.total_exposure_score##i.age_binary c.total_exposure_score##i.region_class2 `controls' || id_num:
estimates store model_full

* Likelihood ratio test
lrtest model_full model_constrained
local p_lrt = r(p)
local chi2_lrt = r(chi2)
local df_lrt = r(df)
di "Likelihood ratio test: χ²(`df_lrt') = `chi2_lrt', p = `p_lrt'"

* 8. Build stratified results dataset -------------------------------------------------------
di as txt "=== 8. Build stratified results dataset ==="

clear
set obs 12

gen stratify_type  = ""
gen stratify_level = ""
gen coef  = .
gen se    = .
gen ci_low  = .
gen ci_high = .
gen p_value = .

* Populate stratification labels
* sex分层 1-3
forval i = 1/3 {
    replace stratify_type = "sex" if _n == `i'
    replace stratify_level = "Full sample" if _n == 1
    replace stratify_level = "Male"   if _n == 2
    replace stratify_level = "Female"   if _n == 3
}

* age分层 4-6
forval i = 4/6 {
    replace stratify_type = "age" if _n == `i'
    replace stratify_level = "Full sample" if _n == 4
    replace stratify_level = "<=65 yrs"  if _n == 5
    replace stratify_level = ">65 yrs"  if _n == 6
}

* Region stratification 7-9
forval i = 7/9 {
    replace stratify_type = "区域" if _n == `i'
    replace stratify_level = "Full sample" if _n == 7
    replace stratify_level = "North"   if _n == 8
    replace stratify_level = "South"   if _n == 9
}

* Populate coefficients and standard errors
* Full sample
replace coef = `score_coef_full'  if stratify_level == "Full sample"
replace se   = `score_se_full'    if stratify_level == "Full sample"

* sex分层
replace coef = `score_coef_male'  if stratify_level == "Male"
replace se   = `score_se_male'    if stratify_level == "Male"

replace coef = `score_coef_female' if stratify_level == "Female"
replace se   = `score_se_female'   if stratify_level == "Female"

* age分层
replace coef = `score_coef_young' if stratify_level == "<=65 yrs"
replace se   = `score_se_young'   if stratify_level == "<=65 yrs"

replace coef = `score_coef_old'   if stratify_level == ">65 yrs"
replace se   = `score_se_old'     if stratify_level == ">65 yrs"

* Region stratification
replace coef = `score_coef_north' if stratify_level == "North"
replace se   = `score_se_north'   if stratify_level == "North"

replace coef = `score_coef_south' if stratify_level == "South"
replace se   = `score_se_south'   if stratify_level == "South"

* Calculate confidence intervals and p-values
replace ci_low  = coef - 1.96 * se
replace ci_high = coef + 1.96 * se
replace p_value = 2 * (1 - normal(abs(coef/se)))

* Add interaction p-values
gen interaction_p = .
replace interaction_p = `p_interaction_gender' if stratify_type == "sex"
replace interaction_p = `p_interaction_age'    if stratify_type == "age"  
replace interaction_p = `p_interaction_region' if stratify_type == "区域"

* 添加Sample size
gen n = .
replace n = `n_full' if stratify_level == "Full sample"
replace n = `n_male' if stratify_level == "Male"
replace n = `n_female' if stratify_level == "Female"
replace n = `n_young' if stratify_level == "<=65 yrs"
replace n = `n_old' if stratify_level == ">65 yrs"
replace n = `n_north' if stratify_level == "North"
replace n = `n_south' if stratify_level == "South"

* Plot position
gen forest_pos = _n

save "output/tables/异质性分析完整结果.dta", replace

* 9. Heterogeneity quantification ---------------------------------------------------------
di as txt "=== 9. Heterogeneity quantification ==="

* Check if metan installed
capture which metan
if _rc {
    di "metan not installed, installing..."
    ssc install metan
}

* Build effect size dataset for meta-analysis
clear
set obs 6

gen subgroup = ""
gen coef = .
gen se = .
gen n = .

replace subgroup = "Male" if _n == 1
replace coef = `score_coef_male' if _n == 1
replace se = `score_se_male' if _n == 1
replace n = `n_male' if _n == 1

replace subgroup = "Female" if _n == 2
replace coef = `score_coef_female' if _n == 2
replace se = `score_se_female' if _n == 2
replace n = `n_female' if _n == 2

replace subgroup = "<=65 yrs" if _n == 3
replace coef = `score_coef_young' if _n == 3
replace se = `score_se_young' if _n == 3
replace n = `n_young' if _n == 3

replace subgroup = ">65 yrs" if _n == 4
replace coef = `score_coef_old' if _n == 4
replace se = `score_se_old' if _n == 4
replace n = `n_old' if _n == 4

replace subgroup = "North" if _n == 5
replace coef = `score_coef_north' if _n == 5
replace se = `score_se_north' if _n == 5
replace n = `n_north' if _n == 5

replace subgroup = "South" if _n == 6
replace coef = `score_coef_south' if _n == 6
replace se = `score_se_south' if _n == 6
replace n = `n_south' if _n == 6

* 计算Heterogeneity statistics
metan coef se, random
local i2 = r(I_sq)
local tau2 = r(tau2)
local p_heterogeneity = r(p_het)
local q_het = r(chi2)

di "Heterogeneity statistics:"
di "I² = `i2'%"
di "τ² = `tau2'"
di "Q统计量 = `q_het'"
di "Heterogeneity test p-value = `p_heterogeneity'"

* Save heterogeneity statistics
gen i2 = `i2'
gen tau2 = `tau2'
gen p_heterogeneity = `p_heterogeneity'
gen q_statistic = `q_het'

save "output/tables/效应量数据集.dta", replace




* 10. Result visualization ------------------------------------------------------------
di as txt "=== 10. Result visualization ==="

* 10.1 Forest plot - fixed version
use "output/tables/异质性分析完整结果.dta", clear

* 只保留分层结果，排除重复的Full sample
keep if stratify_level != "Full sample"

* Drop existing forest_pos, regenerate
capture drop forest_pos
gen forest_pos = _n

* Check data order and labels
list stratify_type stratify_level forest_pos coef ci_low ci_high

* Re-retrieve interaction p-values (local macros may be lost)
sum interaction_p if stratify_type == "sex"
local p_interaction_gender = r(mean)

sum interaction_p if stratify_type == "age"  
local p_interaction_age = r(mean)

sum interaction_p if stratify_type == "区域"
local p_interaction_region = r(mean)

* Format interaction p-values for readability
local p_gender_fmt = string(`p_interaction_gender', "%9.3f")
local p_age_fmt = string(`p_interaction_age', "%9.3f") 
local p_region_fmt = string(`p_interaction_region', "%9.3f")

* Draw forest plot
#delimit ;
twoway 
    (scatter forest_pos coef, 
        mcolor(blue) msymbol(O) msize(medium))
    (rspike ci_high ci_low forest_pos, 
        horizontal lcolor(blue) lwidth(medium)),
    ylabel(
        1 "Male" 
        2 "Female" 
        3 "<=65 yrs" 
        4 ">65 yrs" 
        5 "North" 
        6 "South",
        angle(0) labsize(small) notick)
    ytitle("")
    xtitle("Mixed exposure score coefficient estimate (95% CI)", size(small))
    title("Heterogeneity analysis of EVI and O3 mixed exposure score", size(medium))
    subtitle("按sex、age和Region stratification", size(small))
    xline(0, lpattern(dash) lcolor(red))
    text(0.5 -1.2 "交互作用p值: sex=`p_gender_fmt', age=`p_age_fmt', 区域=`p_region_fmt'", 
         placement(e) size(vsmall) just(left))
    legend(off)
    plotregion(margin(small))
    graphregion(color(white))
    yscale(reverse range(0.5 6.5))
    xlabel(-1.5(0.5)1.5)
    name(heterogeneity_forest, replace);
#delimit cr

graph export "output/figures/heterogeneity_forest.png", replace width(2000) height(1400)
di "Forest plot exported: 异质性分析森林图.png"




* 11. Export results as CSV tables -----------------------------------------------------

di as txt "=== 11. Export results as CSV tables（稳健版本）==="

* 11.1 Reload complete results and generate all necessary variables
use "output/tables/异质性分析完整结果.dta", clear

* Calculate SE (back-calculate from CI)
gen se_calculated = (ci_high - ci_low) / (2 * 1.96)

* Generate significance asterisks
gen significance = ""
replace significance = "***" if p_value < 0.001
replace significance = "**" if p_value >= 0.001 & p_value < 0.01
replace significance = "*" if p_value >= 0.01 & p_value < 0.05
replace significance = "" if p_value >= 0.05

gen interaction_significance = ""
replace interaction_significance = "***" if interaction_p < 0.001
replace interaction_significance = "**" if interaction_p >= 0.001 & interaction_p < 0.01
replace interaction_significance = "*" if interaction_p >= 0.01 & interaction_p < 0.05
replace interaction_significance = "" if interaction_p >= 0.05

* Generate display text
gen coefficient_display = string(coef, "%9.3f") + " (" + string(ci_low, "%9.3f") + ", " + string(ci_high, "%9.3f") + ")" + significance
gen p_value_display = string(p_value, "%9.3f")
gen interaction_p_display = string(interaction_p, "%9.3f") + interaction_significance

* 11.2 Export main results summary table
preserve
keep stratify_type stratify_level n coefficient_display p_value_display interaction_p_display
export delimited using "output/tables/交互项_Stratified analysis_结果汇总.csv", replace
di "Exported: 交互项_Stratified analysis_结果汇总.csv"
restore

* 11.3 Export detailed statistical results
preserve
keep stratify_type stratify_level n coef se_calculated ci_low ci_high p_value interaction_p
export delimited using "output/tables/交互项_Stratified analysis_详细结果.csv", replace
di "Exported: 交互项_Stratified analysis_详细结果.csv"
restore

* 11.4 Create interaction results table (simplified)
di "Creating interaction results table..."

* 重新运行Interaction Analysis来获取p值
preserve
use "DATA/analysis_ready/analysis_dataset.dta", clear
keep if age > 45 & age < 90
gen ln_hh_consumption_total = log(hh_consumption_total)
gen age_binary = 0 if age <= 65
replace age_binary = 1 if age > 65 & age < .
xtile EVI_quartile = EVI_1KM_mean4y, nq(4)
xtile O3_quartile = O3_1KM_mean4y, nq(4)
gen EVI_score = EVI_quartile - 1
gen O3_score = 3 - (O3_quartile - 1)
gen total_exposure_score = EVI_score + O3_score
pca TMP PER SNOW DEM lat
predict climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5, score

local controls "c.age i.sex i.education c.ln_hh_consumption_total i.retirement i.pension_insurance i.smoking_current i.alcohol_current climate_pc1 climate_pc2 climate_pc3 climate_pc4 climate_pc5 c.bmi i.hypertension i.dyslipidemia i.year"

* sex交互作用
mixed BAA_PhenoAge c.total_exposure_score##i.sex `controls' || id_num:
testparm i.sex#c.total_exposure_score
local gender_p = r(p)

* age交互作用
mixed BAA_PhenoAge c.total_exposure_score##i.age_binary `controls' || id_num:
testparm i.age_binary#c.total_exposure_score
local age_p = r(p)

* Region interaction
mixed BAA_PhenoAge c.total_exposure_score##i.region_class2 `controls' || id_num:
testparm i.region_class2#c.total_exposure_score
local region_p = r(p)

restore

* 创建Interaction results表格
clear
set obs 3

gen interaction_type = ""
gen p_value = .
gen interpretation = ""

replace interaction_type = "sex × 混合暴露得分" if _n == 1
replace p_value = `gender_p' if _n == 1
replace interpretation = "检验混合暴露得分效应在男女sex间是否存在差异" if _n == 1

replace interaction_type = "age × 混合暴露得分" if _n == 2  
replace p_value = `age_p' if _n == 2
replace interpretation = "检验混合暴露得分效应在不同age组间是否存在差异" if _n == 2

replace interaction_type = "Region x Mixed Exposure Score" if _n == 3
replace p_value = `region_p' if _n == 3
replace interpretation = "Test whether mixed exposure score effect differs by north-south region" if _n == 3

* Add significance markers
gen significance = ""
replace significance = "***" if p_value < 0.001
replace significance = "**" if p_value >= 0.001 & p_value < 0.01
replace significance = "*" if p_value >= 0.01 & p_value < 0.05
replace significance = "NS" if p_value >= 0.05

export delimited using "output/tables/Interaction Analysis结果.csv", replace
di "Exported: Interaction Analysis结果.csv"

* 11.5 导出Stratified analysis结果
use "output/tables/异质性分析完整结果.dta", clear

* 筛选出分层结果（排除Full sample）
keep if stratify_level != "Full sample"

* Recalculate standard error
gen se_calculated = (ci_high - ci_low) / (2 * 1.96)

* Create English variable names for stratification
gen stratify_var_en = ""
replace stratify_var_en = "Gender" if stratify_type == "sex"
replace stratify_var_en = "Age" if stratify_type == "age" 
replace stratify_var_en = "Region" if stratify_type == "区域"

gen subgroup_en = ""
replace subgroup_en = "Male" if stratify_level == "Male"
replace subgroup_en = "Female" if stratify_level == "Female"
replace subgroup_en = "≤65 years" if stratify_level == "<=65 yrs"
replace subgroup_en = ">65 years" if stratify_level == ">65 yrs" 
replace subgroup_en = "Northern" if stratify_level == "North"
replace subgroup_en = "Southern" if stratify_level == "South"

export delimited stratify_var_en subgroup_en n coef se_calculated ci_low ci_high p_value using "output/tables/Stratified analysis结果.csv", replace
di "Exported: Stratified analysis结果.csv"

* 11.6 Create comprehensive analysis summary report
di "正在Create comprehensive analysis summary report..."

* Get heterogeneity statistics
capture confirm file "output/tables/效应量数据集.dta"
if _rc == 0 {
    use "output/tables/效应量数据集.dta", clear
    sum i2
    local i2_value = string(r(mean), "%9.1f") + "%"
    sum tau2
    local tau2_value = string(r(mean), "%9.3f")
    sum p_heterogeneity
    local p_het_value = string(r(mean), "%9.3f")
}
else {
    local i2_value = "NA"
    local tau2_value = "NA" 
    local p_het_value = "NA"
}

* 获取总Sample size
use "DATA/analysis_ready/analysis_dataset.dta", clear
count if age > 45 & age < 90
local n_total = r(N)

clear
set obs 10

gen category = ""
gen variable = "" 
gen value = ""
gen description = ""

* Interaction results
replace category = "Interaction Analysis" if _n == 1
replace variable = "sex × 暴露得分" if _n == 1
replace value = string(`gender_p', "%9.3f") if _n == 1
replace description = "sex与混合暴露得分的交互作用p值" if _n == 1

replace category = "Interaction Analysis" if _n == 2
replace variable = "age × 暴露得分" if _n == 2  
replace value = string(`age_p', "%9.3f") if _n == 2
replace description = "age与混合暴露得分的交互作用p值" if _n == 2

replace category = "Interaction Analysis" if _n == 3
replace variable = "区域 × 暴露得分" if _n == 3
replace value = string(`region_p', "%9.3f") if _n == 3
replace description = "Region x Mixed Exposure Score interaction p-value" if _n == 3

* Heterogeneity statistics
replace category = "Heterogeneity statistics" if _n == 4
replace variable = "I²" if _n == 4
replace value = "`i2_value'" if _n == 4
replace description = "Degree of heterogeneity" if _n == 4

replace category = "Heterogeneity statistics" if _n == 5
replace variable = "τ²" if _n == 5
replace value = "`tau2_value'" if _n == 5
replace description = "Between-group variance" if _n == 5

replace category = "Heterogeneity statistics" if _n == 6
replace variable = "Heterogeneity test p-value" if _n == 6
replace value = "`p_het_value'" if _n == 6
replace description = "Overall heterogeneity test" if _n == 6

* Sample information
replace category = "Sample information" if _n == 7
replace variable = "总Sample size" if _n == 7
replace value = string(`n_total', "%9.0f") if _n == 7
replace description = "分析总Sample size" if _n == 7

* Analysis information
replace category = "Analysis information" if _n == 8
replace variable = "Analysis time" if _n == 8
replace value = c(current_date) + " " + c(current_time)
replace description = "Analysis execution time" if _n == 8

* Score description
replace category = "Score description" if _n == 9
replace variable = "EVI scoring rule" if _n == 9
replace value = "Q1=0pts, Q2=1pt, Q3=2pts, Q4=3pts"
replace description = "EVI quartile scoring rule" if _n == 9

replace category = "Score description" if _n == 10
replace variable = "O3 scoring rule" if _n == 10
replace value = "Q1=3pts, Q2=2pts, Q3=1pt, Q4=0pts"
replace description = "O3 quartile scoring rule" if _n == 10

export delimited using "output/tables/综合分析汇总报告.csv", replace
di "Exported: 综合分析汇总报告.csv"

* 12. Output file summary
di _n
di as txt "================================================================"
di as txt "              Output file summary"
di as txt "================================================================"
di "生成的CSV表格文件:"
di "1. 交互项_Stratified analysis_结果汇总.csv - 主要结果汇总"
di "2. 交互项_Stratified analysis_详细结果.csv - 详细统计结果"
di "3. Interaction Analysis结果.csv         - 交互作用专用表格"
di "4. Stratified analysis结果.csv             - Stratified analysis专用表格"
di "5. 综合分析汇总报告.csv         - 综合统计报告"
di as txt "================================================================"
di "All tables successfully exported to working directory"
di as txt "================================================================"

di _n
di as txt "=== Result export complete！ ==="