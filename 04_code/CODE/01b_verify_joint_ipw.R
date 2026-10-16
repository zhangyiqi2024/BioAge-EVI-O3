# Exact reconstruction of the corrected archived joint propensity weights.
# Covariate coding is numeric, matching the July 2026 correction.
# Run from 04_code/: Rscript CODE/01b_verify_joint_ipw.R
# Default verifies the delivered input; --write-candidate writes a separate DTA.
suppressPackageStartupMessages({library(haven);library(nnet)})
set.seed(20240101)
d <- as.data.frame(read_dta("DATA/analysis_ready/analysis_dataset.dta"))
b <- d[d$year == 2011,]
stopifnot(!anyDuplicated(b$id_num))
checks <- list()
for(g in c("joint_exp_1KM","joint_exp_500m","joint_exp_1KM_mean4y_base")) {
 w <- sub("joint_exp","ipw_joint",g)
 cols <- c(g,"age","sex","education","urban_rural","marital_status")
 bd <- b[complete.cases(b[,cols]),]
 form <- as.formula(paste("factor(",g,") ~ age + sex + education + urban_rural + marital_status"))
 fit <- multinom(form,data=bd,trace=FALSE,MaxNWts=10000,maxit=500)
 if(fit$convergence != 0) stop("Joint propensity model did not converge: ",g)
 prob <- predict(fit,type="probs")
 selected <- prob[cbind(seq_len(nrow(bd)),match(as.character(bd[[g]]),colnames(prob)))]
 stopifnot(all(is.finite(selected)),all(selected>0),all(selected<=1))
 weights <- 1/selected
 limits <- quantile(weights,c(.01,.99),type=7)
 weights <- pmin(pmax(weights,limits[1]),limits[2])
 aligned <- weights[match(d$id_num,bd$id_num)]
 delta <- max(abs(aligned-d[[w]]),na.rm=TRUE)
 missing_ok <- identical(is.na(aligned),is.na(d[[w]]))
 checks[[g]] <- data.frame(weight=w,baseline_n=nrow(bd),max_absolute_difference=delta,
                          missingness_equal=missing_ok,pass=delta<1e-10 && missing_ok)
 d[[w]] <- aligned
}
dir.create("output/tables",recursive=TRUE,showWarnings=FALSE)
result <- do.call(rbind,checks)
write.csv(result,"output/tables/joint_ipw_reconstruction_verification.csv",row.names=FALSE)
print(result,row.names=FALSE)
if("--write-candidate" %in% commandArgs(trailingOnly=TRUE))
 write_dta(d,"DATA/analysis_ready/analysis_dataset_ipw_reconstructed.dta")
if(!all(result$pass)) stop("Archived joint weight verification failed.")

