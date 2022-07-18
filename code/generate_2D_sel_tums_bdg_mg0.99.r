require(foreach)
source('setup_rscripts_and_sges.r')


subdir_name <- "2D_sel_bdg_highermu"

system(paste0("mkdir ",subdir_name) )
system(paste0("mkdir ../out/",subdir_name) )

dvals <- seq(0, 0.8, by = 0.1)
dval_string <- paste(dvals, collapse = ",")

reps <- 8

foreach(d = dvals)%do%{
  foreach(i = c(1:3))%do%{
	  
	  scriptname <- paste0(subdir_name,'/d',d ,'_ind', i,'.r')

	  #First, write an R script that describes how to write the xml and run the doc
	  rscript <- (generate_Rscript(subdir_name, d, 2, 0.99, 10000, i, 2, 1.1, 1, reps))
	  write(rscript, scriptname)

	  sgescript <- generate_sgescript(24 + d * 10, scriptname)
	  write(sgescript, gsub("\\.r", "\\.sge", scriptname))

  }
}

scriptname <- paste0(subdir_name,'/dmult_ind1.r')	  
rscript <- (generate_Rscript(subdir_name, dval_string, 2, 0.99, 10000, 1, 2, 1.1, 1, 1))
write(rscript, scriptname)
sgescript <- generate_sgescript(36, scriptname)
write(sgescript, gsub("\\.r", "\\.sge", scriptname))


fi <- list.files(subdir_name)
write(paste("#!/usr/bin/sh\n\n", paste0(paste0("qsub ", subdir_name,"/", fi[grep("sge", fi)]), collapse = "\n"), sep = "\n"), paste0("run_", subdir_name, ".sh"))

