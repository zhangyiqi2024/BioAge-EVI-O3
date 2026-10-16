# Process-isolated complete R pipeline. Run from 04_code/.
if(basename(getwd())=="CODE") setwd("..")
if(!file.exists("DATA/analysis_ready/analysis_dataset.dta")) stop("Restricted input missing; see DATA/BUILD_SCOPE.md")
dir.create("logs",recursive=TRUE,showWarnings=FALSE)
dir.create("output/tables",recursive=TRUE,showWarnings=FALSE)
dir.create("output/figures",recursive=TRUE,showWarnings=FALSE)
exe <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
# Dependencies precede figure and benchmark consumers.
numbers <- c(3,4,8,43,45,9:19,27,30,31,34:36,38:42,44,46:54)
all_scripts <- list.files("CODE",pattern="\\.R$")
scripts <- c("00_check_environment.R","01b_verify_joint_ipw.R",
 unlist(lapply(numbers,function(n)all_scripts[startsWith(all_scripts,sprintf("%02d_",n))])))
guards <- list(
 "03_cohort_flow_replication.R"="DATA/restricted/biomarker_merged.dta",
 "41_attrition_analysis.R"=c("DATA/restricted/biomarker_merged.dta","DATA/restricted/phenoage_qc.dta","DATA/restricted/master_merged.dta"),
 "44_outcome_validation_splitsample.R"="DATA/restricted/biomarker_selected.dta",
 "46_selection_exposure_table.R"=c("DATA/restricted/merged_all_waves.dta","DATA/restricted/biomarker_merged.dta"))
history <- data.frame(script=character(),status=character(),elapsed_seconds=numeric())
start <- Sys.time()
for(s in scripts) {
 skip <- s %in% names(guards) && !all(file.exists(guards[[s]]))
 if(s=="38_phenoage_validation.R") skip <-
  !file.exists("DATA/analysis_ready/analysis_dataset_mortality_baseline.dta") &&
  !all(file.exists(c("DATA/restricted/charls/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta","DATA/restricted/charls/2020charls/Exit_Module.dta","DATA/restricted/charls/2020charls/Sample_Infor.dta")))
 if(skip) {
  cat("[SKIP]",s,": restricted inputs absent\n")
  history <- rbind(history,data.frame(script=s,status="SKIP_RESTRICTED",elapsed_seconds=0))
  next
 }
 cat("[RUN]",s,"\n");flush.console()
 before <- Sys.time()
 code <- system2(exe,args=c("--vanilla",shQuote(file.path("CODE",s))),stdout=file.path("logs",paste0(s,".stdout.log")),stderr=file.path("logs",paste0(s,".stderr.log")))
 history <- rbind(history,data.frame(script=s,status=if(code==0)"PASS" else "FAIL",elapsed_seconds=as.numeric(difftime(Sys.time(),before,units="secs"))))
 write.csv(history,"logs/r_pipeline_status.csv",row.names=FALSE)
 if(code!=0) stop("Step failed: ",s,"; see its stderr log")
}
write.csv(history,"logs/r_pipeline_status.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"logs/master_session_info.txt")
cat("R pipeline completed in",as.numeric(difftime(Sys.time(),start,units="mins")),"minutes.\n")

