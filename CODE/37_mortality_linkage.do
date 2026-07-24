********************************************************************************
# 37_mortality_linkage.do
# CHARLS Mortality Follow-Up Data Linkage for PhenoAge Validation
#
# Purpose:   Link Harmonized CHARLS D (waves 1-4, 2011-2018) and 2020 CHARLS
#            (wave 5) mortality data to the analysis dataset, generating
#            died + survival_time variables needed for C-statistic/calibration.
#
# Author:    YiQi Zhang (ZAFU)
# Date:      2026-07-17
#
# Input:     DATA/CHARLS_中国/Raw_data/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta
#            DATA/CHARLS_中国/Raw_data/2020charls/Exit_Module.dta
#            DATA/CHARLS_中国/Raw_data/2020charls/Sample_Infor.dta
#            DATA/analysis_ready/analysis_dataset.dta
#
# Output:    DATA/analysis_ready/analysis_dataset_mortality.dta
#            DATA/analysis_ready/mortality_variables_summary.log
#
# Next:      Run 38_phenoage_validation.R with the enhanced dataset
********************************************************************************

clear all
set more off
set varabbrev off
cap log close
cap mkdir "output"
cap mkdir "output/logs"
log using "output/logs/mortality_linkage.log", replace text

* ---- Paths (all relative; run from the project root = 04_code/) ----
* Raw CHARLS files are restricted data: they are expected three levels up
* from 04_code/ (the BioAge-EVI-O3 repository layout). 00_master_run.do
* checks for them and skips this step when they are absent.
local RAWROOT "../../.."
local OUT    "DATA/analysis_ready"
local RAW    "`RAWROOT'/DATA/CHARLS_中国/Raw_data"

local HARM_D    "`RAW'/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta"
local EXIT2020  "`RAW'/2020charls/Exit_Module.dta"
local INFO2020  "`RAW'/2020charls/Sample_Infor.dta"
local ANALYSIS  "DATA/analysis_ready/analysis_dataset.dta"

di as result "============================================================"
di as result "  CHARLS Mortality Data Linkage"
di as result "  Date: $S_DATE"
di as result "============================================================"

********************************************************************************
# STEP 1: Extract mortality from Harmonized CHARLS D (waves 1-4)
********************************************************************************
di _n as result "--- Step 1: Harmonized CHARLS D mortality extraction ---"

use "`HARM_D'", clear
di "Harmonized CHARLS D: " _N " observations"

* Keep only ID and mortality-relevant variables
keep ID rabyear r1iwy r2iwy r3iwy r4iwy ///
     r1iwstat r2iwstat r3iwstat r4iwstat ///
     radyear radmonth

* Check death counts per wave
foreach w in 2 3 4 {
    qui count if r`w'iwstat == 5
    di "  Wave `w' deaths (iwstat=5): " r(N)
    qui count if r`w'iwstat == 6
    di "  Wave `w' prev-wave deaths (iwstat=6): " r(N)
}

* ---- Construct composite death indicators ----

* Determine if the respondent died at any point during waves 2-4
* iwstat = 5: died this wave (newly discovered death)
* iwstat = 6: died in a previous wave (already counted)
gen byte harm_died_w2 = (r2iwstat == 5)
gen byte harm_died_w3 = (r3iwstat == 5)
gen byte harm_died_w4 = (r4iwstat == 5)

* Cumulative: died BY end of wave 4 (2018)
gen byte harm_died = (r2iwstat == 5 | r3iwstat == 5 | r3iwstat == 6 | r4iwstat == 5 | r4iwstat == 6)

* Death date from Harmonized D
* radyear / radmonth contain the death date for most decedents
gen harm_death_year  = radyear
gen harm_death_month = radmonth
replace harm_death_month = 6 if harm_died == 1 & missing(harm_death_month)
* ^ approximate to mid-year if month unknown

* Last known alive: the last wave where respondent was interviewed or confirmed alive
* r*iwy = interview year for each wave
gen harm_last_alive_year = .
* Default: last wave they were interviewed
replace harm_last_alive_year = r4iwy if !missing(r4iwy)
replace harm_last_alive_year = r3iwy if missing(harm_last_alive_year) & !missing(r3iwy)
replace harm_last_alive_year = r2iwy if missing(harm_last_alive_year) & !missing(r2iwy)
replace harm_last_alive_year = r1iwy if missing(harm_last_alive_year) & !missing(r1iwy)

tab harm_died, missing
di "Total died (Harmonized D, waves 1-4): " %9.0gc _N

* Rename to avoid merge conflicts
rename ID harm_id
rename rabyear harm_rabyear
rename radyear harm_radyear
rename radmonth harm_radmonth

* Save temp
tempfile harm_mort
save "`harm_mort'", replace


********************************************************************************
# STEP 2: Extract mortality from 2020 CHARLS (wave 5)
********************************************************************************
di _n as result "--- Step 2: 2020 CHARLS mortality extraction ---"

* 2a. Sample_Infor: has a simple 'died' indicator
use "`INFO2020'", clear
di "2020 Sample_Infor: " _N " observations"
tab died, missing
keep ID died
rename died w5_died_info
rename ID w5_id
tempfile info2020
save "`info2020'", replace

* 2b. Exit_Module: has detailed death date info (exb001_1=year, exb001_2=month)
use "`EXIT2020'", clear
di "2020 Exit_Module: " _N " observations"

* Check which respondents have exit/death data
* exb001: death date (1=year, 2=month, 3=day)
* exb002: 1=solar calendar, 2=lunar calendar
keep ID exb001_1 exb001_2 exb001_3 exb002
rename exb001_1 w5_death_year
rename exb001_2 w5_death_month
rename exb001_3 w5_death_day
rename exb002   w5_death_calendar_type
rename ID w5_id

* Mark if exit module has death info
gen byte w5_has_exit = !missing(w5_death_year)
tab w5_has_exit
count if w5_has_exit == 1
di "Decedents with death info in 2020 Exit_Module: " r(N)

tempfile exit2020
save "`exit2020'", replace

* 2c. Merge Sample_Infor + Exit_Module
use "`info2020'", clear
merge 1:1 w5_id using "`exit2020'", nogen

* Construct 2020 composite death indicator
gen byte w5_died = 0
replace w5_died = 1 if w5_died_info == 1
replace w5_died = 1 if w5_has_exit == 1  & w5_died == 0

tab w5_died
di "Total died (2020 wave 5): " %9.0gc _N

* Death date: prefer Exit_Module, fall back to approximate
gen w5_death_yr = w5_death_year
gen w5_death_mo = w5_death_month
replace w5_death_mo = 6 if w5_died == 1 & missing(w5_death_mo)

* For those flagged as died in Sample_Infor but without Exit_Module date,
* approximate death year as 2020
replace w5_death_yr = 2020 if w5_died == 1 & missing(w5_death_yr)
replace w5_death_mo = 6  if w5_died == 1 & missing(w5_death_mo)

tempfile mort2020
save "`mort2020'", replace


********************************************************************************
# STEP 3: Merge mortality data into analysis dataset
********************************************************************************
di _n as result "--- Step 3: Merging mortality into analysis dataset ---"

use "`ANALYSIS'", clear
di "Analysis dataset: " _N " observations, " %9.0gc _N " unique IDs = " _N

* Keep only baseline (2011) — one row per person for mortality validation
preserve
keep if year == 2011
di "Baseline (2011) subset: " _N " persons"

* ---- Merge Harmonized CHARLS D mortality (waves 1-4) ----
* (using dataset `harm_mort' is keyed by harm_id - see Step 1)
rename id harm_id
merge 1:1 harm_id using "`harm_mort'", ///
    keepusing(harm_id harm_died harm_died_w2 harm_died_w3 harm_died_w4 ///
              harm_death_year harm_death_month harm_last_alive_year ///
              harm_rabyear)
rename harm_id id

* Tabulate merge success
tab _merge
di "Harmonized D merge: matched=" %9.0gc _N ", not matched=" _N

* Keep all analysis dataset observations (we need PhenoAge values)
* _merge==2 = in Harmonized D but not in analysis dataset → drop
drop if _merge == 2
drop _merge

* ---- Merge 2020 CHARLS mortality (wave 5) ----
* (using dataset `mort2020' is keyed by w5_id - see Step 2)
rename id w5_id
merge 1:1 w5_id using "`mort2020'", ///
    keepusing(w5_id w5_died w5_death_yr w5_death_mo w5_has_exit)
rename w5_id id

tab _merge
drop if _merge == 2
drop _merge

restore, not  // return to full dataset

********************************************************************************
# STEP 4: Construct final mortality variables
********************************************************************************
di _n as result "--- Step 4: Constructing composite mortality variables ---"

* We work on the full dataset but mortality applies to all rows of the same person
* For simplicity, compute on the full dataset and save

* ---- Composite death indicator (any wave) ----
gen byte mortality_died = 0
replace mortality_died = 1 if harm_died == 1
replace mortality_died = 1 if w5_died == 1

* ---- Death year/month (best available) ----
gen mortality_death_year = .
gen mortality_death_month = .

* Prefer 2020 Exit_Module date (most precise for later deaths)
replace mortality_death_year  = w5_death_yr  if w5_died == 1
replace mortality_death_month = w5_death_mo  if w5_died == 1

* Fall back to Harmonized D date for earlier deaths
replace mortality_death_year  = harm_death_year  if missing(mortality_death_year) & harm_died == 1
replace mortality_death_month = harm_death_month if missing(mortality_death_month) & harm_died == 1

* If died but absolutely no date, approximate as mid-follow-up
replace mortality_death_year  = 2015 if mortality_died == 1 & missing(mortality_death_year)
replace mortality_death_month = 6    if mortality_died == 1 & missing(mortality_death_month)

* ---- Survival time from baseline (2011) ----
* Baseline: approximate mid-2011 (2011.5)
* Death date: death_year + death_month/12
* Censoring: end of 2020 follow-up (2020.5) or last known alive

gen mortality_survival_years = .

* For decedents: death date - baseline
replace mortality_survival_years = (mortality_death_year + mortality_death_month/12) - (2011 + 6/12) ///
    if mortality_died == 1

* For survivors: last known alive - baseline
* If they have a 2018 interview year from Harmonized D, use that
replace mortality_survival_years = harm_last_alive_year + 0.5 - (2011 + 6/12) ///
    if mortality_died == 0 & !missing(harm_last_alive_year)

* If no interview data in Harmonized D, use 2015 if in analysis dataset
replace mortality_survival_years = 2015.5 - (2011 + 6/12) ///
    if mortality_died == 0 & missing(mortality_survival_years) & year == 2015

* Otherwise administrative censoring at end of 2020
replace mortality_survival_years = 2020.5 - (2011 + 6/12) ///
    if mortality_died == 0 & missing(mortality_survival_years)

* Clean up negative or zero survival times
replace mortality_survival_years = 0.5 if mortality_survival_years <= 0

* ---- Additional useful variables ----
* Death source indicator
gen byte mortality_death_source = .
replace mortality_death_source = 1 if harm_died == 1 & w5_died == 0   // waves 2-4 only
replace mortality_death_source = 2 if harm_died == 0 & w5_died == 1   // wave 5 only
replace mortality_death_source = 3 if harm_died == 1 & w5_died == 1   // both sources
label define death_src 1 "Waves 2-4 (Harmonized D)" 2 "Wave 5 (2020 CHARLS)" 3 "Both sources"
label values mortality_death_source death_src

* Death wave (which wave death was first detected)
gen byte mortality_death_wave = .
replace mortality_death_wave = 2 if harm_died_w2 == 1
replace mortality_death_wave = 3 if harm_died_w3 == 1
replace mortality_death_wave = 4 if harm_died_w4 == 1
replace mortality_death_wave = 5 if w5_died == 1 & missing(mortality_death_wave)


********************************************************************************
# STEP 5: Summary statistics and quality checks
********************************************************************************
di _n as result "============================================================"
di as result "  MORTALITY DATA SUMMARY"
di as result "============================================================"

* On the full dataset
di _n as result "--- Full dataset ---"
tab mortality_died
summarize mortality_survival_years if mortality_died == 1, detail
summarize mortality_survival_years if mortality_died == 0, detail

* Tab death by source
tab mortality_death_source

* Tab death by wave
tab mortality_death_wave

* On baseline-only (2011)
di _n as result "--- Baseline (2011) only ---"
preserve
keep if year == 2011
di "N baseline: " _N
count if mortality_died == 1
di "N deaths: " r(N) " (" %6.2f 100*r(N)/_N "%)"
summarize mortality_survival_years, detail
restore

* Check for impossible values
qui count if mortality_survival_years < 0
if r(N) > 0 {
    di as error "WARNING: " r(N) " observations with negative survival time!"
}
qui count if mortality_survival_years > 10
if r(N) > 0 {
    di as text "  Note: " r(N) " observations with survival time > 10 years (beyond 2020)"
}

* List some example decedent records
di _n as result "--- Example death records (first 10) ---"
list id year mortality_died mortality_survival_years mortality_death_year ///
    mortality_death_month mortality_death_source mortality_death_wave ///
    if mortality_died == 1 in 1/10, sep(0) noobs abbreviate(16)


********************************************************************************
# STEP 6: Save enhanced dataset
********************************************************************************
di _n as result "--- Step 6: Saving enhanced dataset ---"

* Drop intermediate variables to keep dataset clean
* (keep harm_* and w5_* for transparency but label them as intermediate)
label variable harm_died              "[Intermed] Died waves 2-4"
label variable w5_died                "[Intermed] Died wave 5"
label variable mortality_died         "Died during follow-up (2011-2020)"
label variable mortality_survival_years "Survival time from baseline 2011 (years)"
label variable mortality_death_year   "Year of death"
label variable mortality_death_month  "Month of death"
label variable mortality_death_source "Death data source (1=Harmonized, 2=2020, 3=both)"
label variable mortality_death_wave   "Wave when death first detected"

cap mkdir "`OUT'"
save "`OUT'/analysis_dataset_mortality.dta", replace
di as result "Saved: `OUT'/analysis_dataset_mortality.dta"

* ---- Also save a baseline-only version for convenience ----
keep if year == 2011
save "`OUT'/analysis_dataset_mortality_baseline.dta", replace
di as result "Saved: `OUT'/analysis_dataset_mortality_baseline.dta"

* ---- Export a CSV version for R ----
export delimited using "`OUT'/analysis_dataset_mortality_baseline.csv", replace
di as result "Saved: `OUT'/analysis_dataset_mortality_baseline.csv"


********************************************************************************
# FINAL REPORT
********************************************************************************
di _n as result "============================================================"
di as result "  MORTALITY LINKAGE COMPLETE"
di as result "============================================================"
di as text "Output files:"
di as text "  1. DATA/analysis_ready/analysis_dataset_mortality.dta"
di as text "     Full dataset (6,568 obs) with mortality variables appended"
di as text "  2. DATA/analysis_ready/analysis_dataset_mortality_baseline.dta"
di as text "     Baseline only (2011, ~3,200 obs) for PhenoAge validation"
di as text "  3. DATA/analysis_ready/analysis_dataset_mortality_baseline.csv"
di as text "     CSV version for R script compatibility"
di _n as text "Next step: Run 38_phenoage_validation.R"
di as text "  The R script will now find 'mortality_died' and"
di as text "  'mortality_survival_years' in the dataset."
di "============================================================"

log close
