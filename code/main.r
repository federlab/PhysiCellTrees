
#This code is used to generate the PhysiCell tumors examined in Figure 4 of Lewinsohn et al 2022, doi pending. 

#The PhysiCell section of the code is a minorly modified version of that released by Ghaffarizadeh et al 2018.
#It includes bespoke printing functions, and has mutation capabilities that affect just one daughter cell instead of both
#which required tinkering with some of the core code functionality instead of just creating external classes

#The general structure of this code is as follows:

#The following three scripts create custom directories for the three different simulated models explored
#2D neutral boundary-driven growth
source('generate_2D_neut_tums_bdg.r')
#2D selective boundary-driven growth
source('generate_2D_sel_tums_bdg_mg0.99.r')
#3D neutral boundary-driven growth
source('generate_3D_neutral_tums.r')

#For each of the three models, the above code does three major things:
#1) creates proper output directory structure in PhysiCellTrees/out
#2) creates a directory in PhysiCellTrees/code with an .r and .sge file
#      for running 8 PhysiCell tumors on a cluster node. In total, files
#      are created to run 25 PC tumors at each of 9 death rates 
#      (d = 0, 0.1, 0.2, ... 0.8). Mostly death rates are run together, 
#      but one file "dmult" runs a single iteration of multiple death rates.
#3) creates a file matching run_*.sh with all qsub commands that starts 
#      all the PhysiCell runs

#Note - the 'generate*.r' scripts call/rely on
# > setup_rscripts_and_sges.r
#You'll need to a specify your own directory in the 'generate_R_script' function and set it as 'dirpath'

#At this point the PhysiCell tumors can be generated via the following three commands
system('run_2D_neut_bdg.sh')
system('run_2D_sel_bdg_highermu.sh')
system('run_3D_neut.sh')
#It's possible you may need to manually set the permissions of these files before they will run

#NB - the above code requests a huge amount of cluster time and prints a large amount of data
#This relies on the script
# > 'design_PC_runs'


#Now, the tumors must be sampled according to some sampling scheme, which is specified in the below three scripts
#Ultimately, in the PhysiCell tumors, we only explore one sampling scheme and depth (diversified sampling,
# n = 100), but this would be the place to set up alternative sampling approaches
source('sample_3D_neut.r')
source('sample_2D_neut_bdg.r')
source('sample_2D_sel_bdg_highermu.r')


#These rely on the scripts:
# > sample_tumor.r
# > PC_logs_to_trees.r
