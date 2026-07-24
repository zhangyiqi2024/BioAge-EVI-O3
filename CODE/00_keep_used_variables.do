/*******************************************************************************
Stata calling PowerShell version (No BOM fix)
Automatically extract dataset variables used in A* analysis code and generate a trimmed .dta

Fixes:
  - PowerShell-generated Astar_keep_command.do uses UTF-8 without BOM,
    avoiding Stata error "keep is not a valid command name".
*******************************************************************************/

set seed 20240101

version 18.0
clear all
set more off

global data_dir "DATA"
global data_file "$data_dir\analysis_dataset.dta"

* Directory containing A* code files
global code_dir "."

* PowerShell script path: copy scan_Astar_variables_noBOM.ps1 to data_dir
global ps1_file "$data_dir\scan_Astar_variables_noBOM.ps1"

global out_dir "$data_dir\Astar_variable_screening_outputs"
capture mkdir "$out_dir"

global allvars_txt      "$out_dir\Astar_all_dataset_variables.txt"
global out_dta          "$out_dir\analysis_dataset_Astar_variables_only.dta"
global used_txt         "$out_dir\Astar_used_variables.txt"
global unused_txt       "$out_dir\Astar_unused_variables.txt"
global report_csv       "$out_dir\Astar_variable_source_report.csv"
global code_file_list   "$out_dir\Astar_code_files_scanned.txt"
global keep_cmd         "$out_dir\Astar_keep_command.do"

* Optional: force-keep variables. Modify as needed:
* global force_keep "id id_num year wave"
global force_keep ""

capture confirm file "$ps1_file"
if _rc {
    di as error "Cannot find PowerShell script: $ps1_file"
    di as error "Please copy scan_Astar_variables_noBOM.ps1 to data_dir, or modify global ps1_file."
    exit 601
}

use "$data_file", clear

ds
local allvars `r(varlist)'
local n_all : word count `allvars'
di as text "Original dataset variable count: " as result `n_all'

capture erase "$allvars_txt"
file open fall using "$allvars_txt", write replace text
foreach v of local allvars {
    file write fall "`v'" _n
}
file close fall

* Call PowerShell to scan A* code and generate BOM-free keep command
shell powershell -NoProfile -ExecutionPolicy Bypass -File "$ps1_file" -CodeDir "$code_dir" -OutDir "$out_dir" -ForceKeep "$force_keep"

capture confirm file "$keep_cmd"
if _rc {
    di as error "Keep command file not generated: $keep_cmd"
    exit 602
}

* Count scanned code files and retained variables
local ncode = 0
file open f1 using "$code_file_list", read text
file read f1 line
while r(eof)==0 {
    if strtrim(`"`line'"') != "" local ++ncode
    file read f1 line
}
file close f1

local nused = 0
file open f2 using "$used_txt", read text
file read f2 line
while r(eof)==0 {
    if strtrim(`"`line'"') != "" local ++nused
    file read f2 line
}
file close f2

di as text "Scanned A* code files: " as result `ncode'
di as text "Dataset variables referenced in A* code: " as result `nused'

if `ncode' == 0 {
    di as error "No A*.do / A*.R / A*.r files found in code_dir. Check: $code_dir"
    exit 498
}
if `nused' == 0 {
    di as error "No variables identified. Check whether A* code references current dataset variable names."
    exit 499
}

do "$keep_cmd"

compress
save "$out_dta", replace

di as result "Done."
di as text "Trimmed dataset saved to: " as result "$out_dta"
di as text "A* used variable list: " as result "$used_txt"
di as text "A* unused variable list: " as result "$unused_txt"
di as text "Variable-source report: " as result "$report_csv"
di as text "Scanned code file list: " as result "$code_file_list"
di as text "Auto-generated keep command: " as result "$keep_cmd"

describe
