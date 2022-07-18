

generate_Rscript <- function(dir, deathrate, dimension, mu, cellNumEndCondition, index, driverLim, fitness, threshold, reps = 8){

dimbool <- "\'true\'"
sizeval <- 1000
depthval <- 0
if(dimension == 3){
	     dimbool <- "\'false\'"
	     sizeval <- 400
	     depthval <- sizeval
}

deathrate_dirname <- deathrate
if(grepl(',', deathrate)){ deathrate_dirname <- "mult" }

command <- sprintf('require(foreach)
require(tidyverse)
require(R.matlab)
require(ggtree)
require(phytools)
require(viridis)
require(cowplot)
require(doParallel)

source("PC_logs_to_trees.r")
source("design_PC_runs.r")

dirpath <- "/net/feder/vol1/project/tumor_evolutionary_sims/out/"

daydir <- "%s/d%s_dim%s_w%s_mu%s_N%s_i%s/"

setup_configs_cube(width = %s, height = %s, depth = %s, threads = 1, migration = 0, 
                   fitness = %s, death = c(%s),
                   threshold = %s, mutation_growth = %s, mutation_mig = 1, 
                   time = 200000000, reps = %s, cellNumEndCondition = %s,
                   driverLim = %s, initialCellNumber = 1, daydir, use2D = %s)

runner_restarter(daydir)
process(daydir)
parse_relationships(daydir)', dir, deathrate_dirname, dimension, fitness, mu, cellNumEndCondition, index, sizeval, sizeval, depthval, 
			      fitness, deathrate, threshold, mu, reps, cellNumEndCondition, driverLim, dimbool)

     return(command)
}


generate_sgescript <- function(time, rfile){

return(sprintf('#$ -S /bin/bash
#$ -cwd
#$ -l hostname="sage037|sage038|sage039"
#$ -o slim-test-stdo.output
#$ -e slim-test-stderr.output
#$ -l mfree=1G
#$ -l h_rt=%s:0:0
#$ -pe serial 8

module load modules modules-init modules-gs pcre2/10.39 hdf5/1.10.1 R/4.1.2

Rscript %s', time, rfile))

}


