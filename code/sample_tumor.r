#Figure out growth rate differences over time

require(foreach)
require(tidyverse)
require(concaveman)

source('diversified_sampling.r')

setUpSampleDir <- function(dirToProcess, outName){

    figdir <- paste0("../out/",dirToProcess,"/", outName, "/figs/")

    system(paste0("mkdir ", paste0("../out/",dirToProcess, "/", outName, "")))
    system(paste0("mkdir ", paste0("../out/", dirToProcess, "/", outName, "/growth_rate_differences")))
    system(paste0("mkdir ", paste0("../out/", dirToProcess, "/", outName, "/to_beast_format")))
    system(paste0("mkdir ", paste0("../out/", dirToProcess, "/", outName, "/with_hull_info")))
    system(paste0("mkdir ", figdir))

}

postProcess_sample <- function(dirToProcess, outName){

    growth_rate_diffs <- foreach(fi = list.files(paste0("../out/", dirToProcess, "/", outName, "/growth_rate_differences/")), .combine = "rbind")%do%{

        tbl_df(read.table(
            file = paste0("../out/", dirToProcess, "/", outName, "/growth_rate_differences/", fi), header = TRUE, sep = "\t"))

    }

    pdf(paste0(figdir, "/growth_rate_differences_reps.pdf"), height = 5, width = 5)
    print(growth_rate_diffs %>% separate(base, into = c("type", "m", "w", "d", "t", "mg", "mm", "l", "i", "s"), sep = "_") %>% 
    mutate(d = as.numeric(gsub("d", "", d))) %>% select(file, edgeAssociated, avGrowth, d, i, s) %>% spread(edgeAssociated, avGrowth) %>% 
        rename(center = `FALSE`, edge = `TRUE`) %>% mutate(diff = edge - center) %>% mutate(isna = is.na(diff), diff = ifelse(is.na(diff), 0, diff)) %>%
            filter(isna == FALSE, file < 110) %>% 
            ggplot() + geom_point(aes(x = file, y = diff/edge, col = factor(d)), alpha = 0.5) +
                geom_path(aes(x = file, y = diff/edge, col = factor(d), group = paste0(i, s)), alpha = 0.2) + labs(y =  "Birth rate difference (b_edge - b_center)/b_edge", x = "Time", col = "Death rate"))
    dev.off()

    pdf(paste0(figdir, "/growth_rate_differences.pdf"), height = 5, width = 5)
    print(
        growth_rate_diffs %>% separate(base, into = c("type", "m", "w", "d", "t", "mg", "mm", "l", "i", "s"), sep = "_") %>% 
            mutate(d = as.numeric(gsub("d", "", d))) %>% select(file, edgeAssociated, avGrowth, d, i, s) %>% spread(edgeAssociated, avGrowth) %>% 
        rename(center = `FALSE`, edge = `TRUE`) %>% mutate(diff = edge - center) %>% mutate(isna = is.na(diff), diff = ifelse(is.na(diff), 0, diff)) %>%
            filter(isna == FALSE, file < 110)  %>% group_by(s, d) %>% filter(file > quantile(file, 0.1)) %>% summarize(diff = mean(diff/edge)) %>% ggplot() + geom_jitter(aes(x = d, y = diff, col = factor(d)), width = 0.01, height = 0) + labs(x = "Death rate", y = "Birth rate difference (b_edge - b_center)/b_edge", col = "death")
    )
    dev.off()

}




sampleTrees <- function(dirToProcess, sampleSize, diversifiedOrRandom, dimension, outName, dirstub, mu = 1, printFigs = FALSE, printHull = FALSE){

#    foreach(dirstub = list_of_dirs)%do%{

    print(paste0("Processing ", dirstub))

    reldir <- sort(unique(gsub("(_par_child|_child_list)?\\.csv", "", list.files(paste0("../out/",dirToProcess, "/", dirstub, "/processed/")))))

    foreach(base = reldir)%do%{

        print(paste0(" -> ", base))
        
        thres <- as.numeric(gsub(".*_t|_mg.*", "", base))

#Three different input types:
        inf_pos_pres <- tbl_df(read.table(paste0("../out/",dirToProcess, "/", dirstub, "/processed/", base, ".csv"),header = TRUE, sep = ","))

#relationships
        inf_par_child <- tbl_df(read.table(paste0("../out/",dirToProcess, "/", dirstub, "/processed/", base, "_par_child.csv"),header = TRUE, sep = ","))

        inf_pos_pres <- inf_pos_pres %>% group_by(file) %>% mutate(meanx = mean(x), meany = mean(y), meanz = mean(z)) %>% ungroup() %>% 
            mutate(distFromCenter = sqrt((x - meanx)^2 + (y - meany)^2 + (z - meanz)^2))

        mf <- (inf_pos_pres %>% summarize(maxfile = max(file)))$maxfile

        if(diversifiedOrRandom == "diversified"){
            if(dimension == 3){
                sampleIDs <- (sample_physicell_diversified_3D(inf_pos_pres %>% filter(file == mf) %>% mutate(trial = base), base, sampleSize))$ID
            }
            if(dimension == 2){
                sampleIDs <- (sample_physicell_diversified(inf_pos_pres %>% filter(file == mf) %>% mutate(trial = base), base, sampleSize))$ID
            }
        }else{
            sampleIDs <- (inf_pos_pres %>% filter(file == mf))$ID
        }

	#Previously, I was just importing the full relationship tree that I have precomputed
	#Due to memory issues, it's actually easier to just recompute this

       #Let's quickly traverse this list backwards from all our extant lineages
        relInds <- c()

        foreach(id = sampleIDs)%do%{

            while(!is.element(id, relInds) ){

                relInds <- append(relInds, id)
                id <- (inf_par_child %>%
                           filter(index == id))$parent_index

                if(length(id) == 0){ break }
            }
        }

        relInds <- sort(relInds, decreasing = TRUE)

        child_list <- sapply(paste(relInds),function(x) NULL)
       

        for(ind in relInds){

                                        #All nodes are their own descendants
            mychildren <- (inf_par_child %>% filter(parent_index == ind))$index
            relChildren <- mychildren[is.element(mychildren, relInds)]

            if(length(relChildren) == 0) {
                child_list[[paste(ind)]] <- c(ind)
            }else{
                child_list[[paste(ind)]] <- c(as.vector(unlist(child_list[paste(relChildren)])), ind)
            }

        }


	relationshipTree <- list()
        relationshipTree$array_format <- inf_par_child
        relationshipTree$child_list  <- child_list

        treeDat <- sampleTree(relationshipTree, sampleSize, mu, rseed = 5, sampleIDs = sampleIDs)
        tiplabs <- treeDat$muttree$tip.label
 
        labeled_tips <- (inf_pos_pres %>% filter(is.element(ID, as.numeric(tiplabs))) %>% filter(file == mf) %>% 
                     dplyr::select(ID, distFromCenter) %>% rename(label = ID) %>% mutate(label = paste(label)))
       
        if(dimension == 2){

            distToHull <- function(point.x, point.y, point.z, hull){

                minDist = hull %>% mutate(distToPoint = sqrt((x - point.x)^2 + (y - point.y)^2 + (z - point.z)^2)) %>% 
                    summarize(mD = min(distToPoint))

                return(minDist$mD)
               
            }

            inf_pos_pres_dth <- foreach(focal_file = unique(inf_pos_pres$file), .combine = "rbind")%do%{

                tmp <- (inf_pos_pres %>% filter(file == focal_file))
                hull <- tbl_df(concaveman(as.matrix(cbind(tmp$x, tmp$y, tmp$z)))) %>% rename(x = V1, y = V2, z = V3)
                tmp <- tmp %>% rowwise() %>% mutate(distToHull = distToHull(x, y, z, hull))

                return(tmp)
                
            }
    
        }
       
        if(dimension == 3){

            slice_thickness <- 10
	    #            zslices <- unique((endStates %>% filter(trial == base))$zslice)
	    zslices <- sort(unique(floor((inf_pos_pres %>% filter(file == mf))$z /  (slice_thickness * 2) ) * slice_thickness * 2))
            sampleable_cells <- purrr::map(zslices, function(zslice){ inf_pos_pres %>% 
                filter(between(z, zslice - slice_thickness, zslice + slice_thickness)) %>% 
                 mutate(zslice = zslice) } ) %>% 
            bind_rows
          
            distToHull <- function(point.x, point.y, point.z, hull){

                minDist = hull %>% mutate(distToPoint = sqrt((x - point.x)^2 + (y - point.y)^2) + (z - point.z)^2) %>% 
                    summarize(mD = min(distToPoint))

                return(minDist$mD)
                
            }

            inf_pos_pres_dth <- foreach(focal_file = unique(inf_pos_pres$file), .combine = "rbind")%do%{

                zoff <- 20
                tmp <- (inf_pos_pres %>% filter(file == focal_file)) %>% 
                    mutate(zslice = floor(z/zoff)*zoff)

                cavefunc <- function(tmp){

                    tbl_df(concaveman(as.matrix(cbind(tmp$x, tmp$y, tmp$z)))) %>% 
                        rename(x = V1, y = V2, z = V3)
                }

                hull <- tmp %>% group_by(zslice) %>% 
                    do(cavefunc(.)) %>% ungroup()

                tmp <- tmp %>% rowwise() %>% mutate(distToHull = distToHull(x, y, z, hull))

                return(tmp)
                
            }           
        }

        relIndices <- (treeDat$maya_format)$index
        baseGrowth <- min(unique(inf_pos_pres_dth$growth))
        w <- as.numeric(gsub(".*_w|_d.*", "", base))

        suppInfForMayaFormat <- (inf_pos_pres_dth %>% ungroup() %>%
                                     filter(is.element(ID, relIndices)) %>% group_by(ID) %>% #filter(file == min(file))  %>% 
				     filter(file == max(file))  %>%
                                     mutate(numDrivers = round(log(growth/baseGrowth)/log(w))) %>%
                                     select(ID, distFromCenter, distToHull, numDrivers) ) %>% ungroup() %>% rename(index = ID)

        #So, we ideally would smush down this tree
        #We really only one nodes that have two kids and extant samples
        viable_parents <- unique((treeDat$maya_format %>% group_by(parent_index) %>% filter(n() == 2))$parent_index)
        viable_children <- c((treeDat$maya_format %>% filter(extant == TRUE))$index, viable_parents[viable_parents != min(viable_parents)])

        collapsedLineages <- foreach(ind = viable_children, .combine = "rbind")%do%{

            myPar <- (treeDat$maya_format %>% filter(index == ind) )$parent_index
            while(!is.element(myPar, viable_parents)){
                myPar <- (treeDat$maya_format %>% filter(index == myPar) )$parent_index
            }
            treeDat$maya_format %>% filter(index == ind) %>% mutate(parent_index = myPar)
             
        }

        collapsedLineages <- collapsedLineages %>% add_row(treeDat$maya_format %>% filter(index == min(viable_parents)))
        maya_toPrint <- left_join(collapsedLineages %>% select(-X) , suppInfForMayaFormat) %>% arrange(index) 

        finalTP <- max(inf_pos_pres_dth$file) * 360
        maya_toPrint <- maya_toPrint %>% mutate(deathdate = ifelse(extant, finalTP, deathdate))

	#At this point, the locations in maya_toPrint are based on where cells are BORN
	#We need to edit this so that we're recording where cells were AT THE FINAL SAMPLING TIMEPOINT
	
	print(head(inf_pos_pres_dth))
	
	print(dim(maya_toPrint %>% filter(extant == TRUE)))

	print(maya_toPrint %>% filter(extant == TRUE))

	newLabs <- left_join(
                maya_toPrint %>% filter(extant == TRUE) %>% select(-locx, -locy), 
                inf_pos_pres_dth %>% rename(index = ID, newDistToHull = distToHull) %>% ungroup() %>% 
                    filter(file == mf) %>% select(index, newDistToHull, x, y) %>% 
                        rename(locx = x, locy = y) 
            ) 

	maya_toPrint <- bind_rows(maya_toPrint %>% filter(extant == FALSE), 
                                   newLabs) 

        write.table(maya_toPrint,
            file = paste0("../out/", dirToProcess, "/", outName, "/to_beast_format/", base, "_", diversifiedOrRandom, "_m", mu, "_n", sampleSize, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

        avGrowth <- inf_pos_pres_dth %>% ungroup() %>% mutate(edgeAssociated = distToHull <= 20)  %>% 
            mutate(growthRate = growth * (pressure <= thres)) %>%
                group_by(file, edgeAssociated) %>% summarize(avGrowth = mean(growthRate)) %>%
                ggplot() + geom_line(aes(x = file, y = avGrowth, col = edgeAssociated)) +
                labs(x = "Time", y = "Average growth rate")

        write.table(
            inf_pos_pres_dth %>% ungroup() %>% mutate(edgeAssociated = distToHull <= 20)  %>% 
                mutate(growthRate = growth * (pressure <= thres)) %>%
                group_by(file, edgeAssociated) %>% summarize(avGrowth = mean(growthRate), n = n()) %>% mutate(base = base),
            file = paste0("../out/", dirToProcess, "/", outName, "/growth_rate_differences/", base, "_", diversifiedOrRandom, "_m", mu, "_n", sampleSize, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

	if(printHull == TRUE){
	    write.table( inf_pos_pres_dth %>% mutate(base = base),
	              file = paste0("../out/", dirToProcess, "/", outName, "/with_hull_info/", base, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
		      }

        if(printFigs == TRUE){

            g1 <- left_join(treeDat$muttree, labeled_tips) %>% ggtree() + geom_tippoint(aes(col = distFromCenter)) + scale_color_viridis()

            if(dimension == 2){

                g4 <- inf_pos_pres %>% filter(file == mf) %>% 
                    ggplot() + geom_point(aes(x = x, y = y, col = pressure >= thres), size = 1) + 
                    theme_classic() + theme(legend.position = "none") 
                g3 <- inf_pos_pres %>% filter(file == mf)  %>% 
                    ggplot() + geom_point(aes(x = x, y = y, col = factor(growth)), size = 1) +
                    geom_point(aes(x = x, y = y, fill = distFromCenter), col = "red", pch = 21, size = 2, data = 
                               inf_pos_pres %>% filter(is.element(ID, (treeDat$toBD %>% filter(fixed == TRUE))$id) )  %>% 
                                   filter(file == mf)) + theme_classic() +
                                   theme(legend.position = "none") + scale_fill_viridis() + scale_color_grey()
                    
                edge <- inf_pos_pres_dth %>% ggplot() + geom_point(aes(x = x, y = y, col = distToHull <= 20), size = 0.1) + facet_wrap(~file)
                pres <- inf_pos_pres_dth %>% ggplot() + geom_point(aes(x = x, y = y, col = pressure < thres), size = 0.1) + facet_wrap(~file)
                drivs <- inf_pos_pres_dth %>% ggplot() + geom_point(aes(x = x, y = y, col = growth), size = 0.1) + facet_wrap(~file)

                treePlot <- plot_grid( plot_grid(g3, g4, nrow = 2), g1, nrow = 1)
            }

            if(dimension == 3){

              
                g4 <- sampleable_cells %>% filter(file == mf) %>% 
                    ggplot() + geom_point(aes(x = x, y = y, col = pressure >= thres), size = 1) + 
                    theme_classic() + theme(legend.position = "none") + facet_wrap(~zslice, ncol = 1)

                g3 <- sampleable_cells %>% filter(file == mf)  %>% 
                    ggplot() + geom_point(aes(x = x, y = y, col = factor(growth)), size = 1) +
                    geom_point(aes(x = x, y = y, fill = distFromCenter), col = "red", pch = 21, size = 2, data = 
                    sampleable_cells %>% filter(is.element(ID, (treeDat$toBD %>% filter(fixed == TRUE))$id) )  %>% 
                    filter(file == mf)) + theme_classic() +
                    theme(legend.position = "none") + scale_fill_viridis() + 
                    scale_color_grey()+ facet_wrap(~zslice, ncol = 1)

                tps <- round(seq(2, mf, length.out = 5))
                toPlot <- left_join(sampleable_cells %>% filter(is.element(file, tps)),  inf_pos_pres_dth %>% select(-zslice) )

                edge <- toPlot %>%
                    ggplot() + geom_point(aes(x = x, y = y, col = distToHull <= 20), size = 0.1) + facet_grid(zslice ~ file)
                pres <- toPlot %>%
                    ggplot() + geom_point(aes(x = x, y = y, col = pressure < thres), size = 0.1) + facet_grid(zslice ~ file)
                drivs <- toPlot %>%
                    ggplot() + geom_point(aes(x = x, y = y, col = growth), size = 0.1) + facet_grid(zslice ~ file)


                treePlot <- plot_grid( g3, g4, g1, nrow = 1, rel_widths = c(1, 1, 1.5))

            }
           
            pdf(paste0(figdir, "ancestral_locations_",base, "_", diversifiedOrRandom, "_m", mu, "_n", sampleSize, ".pdf"), height = 5, width = 6.5 )
            print(plotGuesses(maya_toPrint))
            dev.off()
            
            pdf(paste0(figdir, "avGrowth_",base,"_", diversifiedOrRandom, "_m", mu, "_n", sampleSize,  ".pdf"), height = 5, width = 6 )
            print(avGrowth)
            dev.off()

            pdf(paste0(figdir, "inferredEdge_",base,"_", diversifiedOrRandom, "_m", mu, "_n", sampleSize,  ".pdf"), height = 8, width = 10 )
            print(edge)
            dev.off()

            pdf(paste0(figdir, "truePressure_",base, "_", diversifiedOrRandom, "_m", mu, "_n", sampleSize, ".pdf"), height = 8, width = 10 )
            print(pres)
            dev.off()

            pdf(paste0(figdir, "drivers_",base,"_", diversifiedOrRandom, "_m", mu, "_n", sampleSize,  ".pdf"), height = 8, width = 10 )
            print(drivs)
            dev.off()

            pdf(paste0(figdir, "trees_",base,"_", diversifiedOrRandom, "_m", mu, "_n", sampleSize,  ".pdf"), height = 8, width = 8 )
            print(treePlot )
            dev.off()
        }
        
    }

}



generate_sampling_rscript <- function(dirToProcess, sampleSize, diversifiedOrRandom, dimension, outName, subdir, mu, printHull = FALSE){

return(sprintf('source("sample_tumor.r")
source("PC_logs_to_trees.r")

sampleTrees("%s", %s, "%s", %s, "%s", "%s", %s, printHull = %s)', 
dirToProcess, sampleSize, diversifiedOrRandom, dimension, outName, subdir, mu, printHull))

}

generate_sampling_sgescript <- function(time, rfile){

return(sprintf('#$ -S /bin/bash
#$ -cwd
#$ -o sampling-stdo.output
#$ -e sampling-stderr.output
#$ -l mfree=4G
#$ -l h_rt=%s:0:0

module load modules modules-init modules-gs pcre2/10.39 hdf5/1.10.1 R/4.1.2

Rscript %s', time, rfile))

}

generate_sample_scheme <- function(dirToProcess, sampleSize, diversifiedOrRandom, dimension, outName, mu, printHull = FALSE){

#first - setup output dir
setUpSampleDir(dirToProcess, outName)

#second - setup dir to hold Rscripts/SGEs
system(paste0("mkdir ", dirToProcess, "/", outName))

#Get a list of all the subdirs to run. We'll need to make an R and sge for each one
subDirs <- list.files(paste0("../out/", dirToProcess))
subDirs <- subDirs[grep("d.*_dim[2-3]*_.*", subDirs)]

foreach(subDir = subDirs)%do%{
	       
	       print(paste0("generating a file for ", subDir))
	       rscript <- generate_sampling_rscript(dirToProcess, sampleSize, diversifiedOrRandom, dimension, outName, subDir, mu, printHull)
	       scriptname <- paste0(dirToProcess,'/',outName,'/', subDir, '.r')
	       sgescript <- generate_sampling_sgescript(24, scriptname)

	       write(rscript, scriptname)
               write(sgescript, gsub("\\.r", "\\.sge", scriptname))

}


subdir_name <- paste0(dirToProcess, "/", outName)
fi <- list.files(subdir_name)
write(paste("#!/usr/bin/sh\n\n", paste0(paste0("qsub ", subdir_name,"/", fi[grep("sge", fi)]), collapse = "\n"), sep = "\n"), paste0("run_", dirToProcess, "_", outName, ".sh"))



}

generate_sample_scheme_varSweep <- function(dirToProcess, sampleSize, diversifiedOrRandom, dimension, outName, mu, printHull = FALSE, 
				             subdirRegExp = NA){

#first - setup output dir
setUpSampleDir(dirToProcess, outName)

#second - setup dir to hold Rscripts/SGEs
system(paste0("mkdir ", dirToProcess, "/", outName))

#Get a list of all the subdirs to run. We'll need to make an R and sge for each one
subDirs <- list.files(paste0("../out/", dirToProcess))
subDirs <- subDirs[grep("d.*_dim[2-3]*_.*", subDirs)]

#Only process directories that match this regular expression 
#for example, to look at a single death rate
if(!is.na(subdirRegExp )){
 subDirs <- subDirs[grepl(subdirRegExp, subDirs) ]
}

foreach(subDir = subDirs)%do%{
       foreach(mu_val = mu)%do%{
	       foreach(sampleSize_val = sampleSize)%do%{
	       
	       print(paste0("generating a file for ", subDir))
	       rscript <- generate_sampling_rscript(dirToProcess, sampleSize_val, diversifiedOrRandom, dimension, outName, subDir, mu_val, printHull)
	       scriptname <- paste0(dirToProcess,'/',outName,'/', subDir, '_n', sampleSize_val,'_mun',mu_val,'.r')
	       sgescript <- generate_sampling_sgescript(24, scriptname)

	       write(rscript, scriptname)
               write(sgescript, gsub("\\.r", "\\.sge", scriptname))
	}
	}
}


subdir_name <- paste0(dirToProcess, "/", outName)
fi <- list.files(subdir_name)
write(paste("#!/usr/bin/sh\n\n", paste0(paste0("qsub ", subdir_name,"/", fi[grep("sge", fi)]), collapse = "\n"), sep = "\n"), paste0("run_", dirToProcess, "_", outName, ".sh"))

}
