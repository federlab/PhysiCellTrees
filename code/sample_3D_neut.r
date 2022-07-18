source('sample_tumor.r')
source('PC_logs_to_trees.r')

dirToProcess <- "3D_neut"
dimension <- 3
sampleSize <- 100
diversifiedOrRandom <- "diversified"
sampleSchemeName <- "diversified_100"
mu <- 1
printHull = TRUE

generate_sample_scheme(dirToProcess, sampleSize, diversifiedOrRandom, dimension, sampleSchemeName, mu, printHull = printHull)



