version 18.0
clear all
set more off
* Run after PYTHON/prepare_analysis_data.py has created the English input.
* This selects fields from an authorised analysis-level source. It does not
* reconstruct the supplied biological-age outcome or spatial exposures.
capture confirm file "DATA/restricted/analysis_input_english.dta"
if _rc {
    display as error "Missing DATA/restricted/analysis_input_english.dta."
    display as error "See DATA/BUILD_SCOPE.md. Do not substitute another outcome."
    exit 601
}
use "DATA/restricted/analysis_input_english.dta", clear
local wanted id year city BA_PhenoAge BAA_PhenoAge wave province_code city_code ///
    sex age education_raw education hukou_type han_ethnicity marital_status urban_rural ///
    employed retirement hypertension dyslipidemia diabetes cancer lung_disease heart_disease ///
    stroke memory_disease arthritis asthma outpatient_visits inpatient_visits pa_vigorous pa_moderate ///
    pa_light days_vigorous days_moderate days_light hobbies self_rated_health adl_difficulty smoking_current ///
    alcohol_current hospital pension_insurance medical_insurance cognitive_ability social_participation hh_consumption_total hh_consumption_percapita ///
    hh_income_total mental_health life_satisfaction pa_moderate_vigorous_weekly metabolic_activity social_isolation1 social_isolation2 disability ///
    bmi outpatient_cost_total bl_tg bl_hdl bl_ldl bl_hbalc region_class1 hu_line ///
    region_class2 region lat PM25_1KM PM10_1KM O3_1KM EVI_1KM PM25_500m ///
    PM10_500m O3_500m EVI_500m WIND SNOW TMP PM25 PM10 ///
    PER O3 EVI DEM PM25_1KM_lag1 PM10_1KM_lag1 O3_1KM_lag1 EVI_1KM_lag1 ///
    PM25_500m_lag1 PM10_500m_lag1 O3_500m_lag1 EVI_500m_lag1 PM25_lag1 PM10_lag1 O3_lag1 EVI_lag1 ///
    O3_1KM_lag2 EVI_1KM_lag2 O3_500m_lag2 EVI_500m_lag2 O3_lag2 EVI_lag2 O3_1KM_lag3 EVI_1KM_lag3 ///
    O3_500m_lag3 EVI_500m_lag3 O3_lag3 EVI_lag3 O3_1KM_lag4 EVI_1KM_lag4 O3_500m_lag4 EVI_500m_lag4 ///
    O3_lag4 EVI_lag4 inversion_days1_lag5 inversion_days2_lag5 steep_slope_pct_25deg_1km_lag5 steep_slope_pct_25deg_500m_lag5 steep_slope_pct_25deg_lag5 steep_slope_pct_20deg_1km_lag5 ///
    steep_slope_pct_20deg_500m_lag5 steep_slope_pct_20deg_lag5 mean_slope_1km_lag5 mean_slope_500m_lag5 mean_slope_lag5 O3_1KM_lag5 EVI_1KM_lag5 O3_500m_lag5 ///
    EVI_500m_lag5 O3_lag5 EVI_lag5 delta_BAA_PhenoAge social_participation_score resilience_A resilience_B days_light_moderate ///
    sleep_quality days_moderate_vigorous days_any_pa pa_light_moderate_weekly pa_any_weekly traj15_PM25_1KM traj15_PM10_1KM traj15_O3_1KM ///
    traj15_NDVI_1KM traj15_FVC_1KM traj15_EVI_1KM traj15_PM25_500m traj15_PM10_500m traj15_O3_500m traj15_NDVI_500m traj15_FVC_500m ///
    traj15_EVI_500m traj15_PM25 traj15_PM10 traj15_O3 traj15_NDVI traj15_FVC traj15_EVI BA_PhenoAge_base_2011 ///
    BAA_PhenoAge_base_2011 EVI_mean2y EVI_mean3y EVI_mean4y EVI_mean5y EVI_500m_mean2y EVI_500m_mean3y EVI_500m_mean4y ///
    EVI_500m_mean5y EVI_1KM_mean2y EVI_1KM_mean3y EVI_1KM_mean4y EVI_1KM_mean5y O3_mean2y O3_mean3y O3_mean4y ///
    O3_mean5y O3_500m_mean2y O3_500m_mean3y O3_500m_mean4y O3_500m_mean5y O3_1KM_mean2y O3_1KM_mean3y O3_1KM_mean4y ///
    O3_1KM_mean5y O3_1KM_mean4y_cat EVI_1KM_mean4y_cat ba_acceleration_A id_num ipw_EVI_500m ipw_EVI_1KM ipw_O3_500m ///
    ipw_O3_1KM ipw_O3_1KM_mean4y_cat ipw_EVI_1KM_mean4y_cat joint_exp_1KM joint_exp_500m joint_exp_1KM_mean4y_base joint_exp_1KM_mean4y_current dataset_type ///
    imputation_date imputation_method ipw_joint_1KM ipw_joint_500m ipw_joint_1KM_mean4y_base

local retained
foreach v of local wanted {
    capture confirm variable `v'
    if !_rc local retained "`retained' `v'"
}
foreach v in id age sex year BAA_PhenoAge BA_PhenoAge {
    confirm variable `v'
}
keep `retained'
capture mkdir "DATA/analysis_ready"
save "DATA/analysis_ready/selected_input.dta", replace
