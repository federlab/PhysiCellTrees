#Figure out growth rate differences over time

require(foreach)
require(tidyverse)
require(concaveman)

#list_of_dirs <- list.files("../out/")[ grep("2022_03_2[6|9]_[0-9]", list.files("../out/"))]
#

list_of_dirs <- c("2022_03_22_01", "2022_03_22_02", "2022_03_23_01", "2022_03_24_01", "2022_03_18_03", "2022_03_18_02", "2022_03_30_01", "2022_03_30_02", "2022_03_30_03")

source('diversified_sampling.r')
source('extra_functions.r')
source('PC_logs_to_trees.r')

#outName <- "noboundarygrowth_withselection_2D"
outName <- "sampling_3D"
system(paste0("mkdir ", paste0("../out/", outName, "")))
system(paste0("mkdir ", paste0("../out/", outName, "/growth_rate_differences")))
system(paste0("mkdir ", paste0("../out/", outName, "/to_beast_format")))

dimension <- 3

foreach(dirstub = list_of_dirs[-1])%do%{

    print(paste0("Processing ", dirstub))

    reldir <- sort(unique(gsub("(_par_child|_child_list)?\\.csv", "", list.files(paste0("../out/", dirstub, "/processed/")))))

    foreach(base = reldir[grep("^samp", reldir)])%do%{

        print(paste0(" -> ", base))
        
        thres <- as.numeric(gsub(".*_t|_mg.*", "", base))

#Three different input types:
        inf_pos_pres <- tbl_df(read.table(paste0("../out/", dirstub, "/processed/", base, ".csv"),header = TRUE, sep = ","))

#relationships
        inf_par_child <- tbl_df(read.table(paste0("../out/", dirstub, "/processed/", base, "_par_child.csv"),header = TRUE, sep = ","))

#descendant
        inf_child_list <- tbl_df(read.table(paste0("../out/", dirstub, "/processed/", base, "_child_list.csv"),header = TRUE, sep = ",")) %>% 
            rename(par = i, desc = x) %>% select(par, desc)

        child_list <- by(inf_child_list, inf_child_list$par, function(x){return(as.numeric(x$desc))})

        relationshipTree <- list()
        relationshipTree$array_format <- inf_par_child
        relationshipTree$child_list  <- child_list

        inf_pos_pres <- inf_pos_pres %>% group_by(file) %>% mutate(meanx = mean(x), meany = mean(y), meanz = mean(z)) %>% ungroup() %>% 
            mutate(distFromCenter = sqrt((x - meanx)^2 + (y - meany)^2 + (z - meanz)^2))

        mf <- (inf_pos_pres %>% summarize(maxfile = max(file)))$maxfile

        if(dimension == 3){

            sampleIDs <- (sample_physicell_diversified_3D(inf_pos_pres %>% filter(file == mf) %>% mutate(trial = base), base, 100))$ID
        }
        if(dimension == 2){
            sampleIDs <- (sample_physicell_diversified(inf_pos_pres %>% filter(file == mf) %>% mutate(trial = base), base, 100))$ID
        }

        treeDat <- sampleTree(relationshipTree, 100, 1, rseed = 5, sampleIDs = sampleIDs)
        tiplabs <- treeDat$muttree$tip.label
 
        labeled_tips <- (inf_pos_pres %>% filter(is.element(ID, as.numeric(tiplabs))) %>% filter(file == mf) %>% 
                     dplyr::select(ID, distFromCenter) %>% rename(label = ID) %>% mutate(label = paste(label)))

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
            
            edge <- inf_pos_pres_dth %>% ggplot() + geom_point(aes(x = x, y = y, col = distToHull <= 20), size = 0.1) + facet_wrap(~file)
            pres <- inf_pos_pres_dth %>% ggplot() + geom_point(aes(x = x, y = y, col = pressure < thres), size = 0.1) + facet_wrap(~file)
            drivs <- inf_pos_pres_dth %>% ggplot() + geom_point(aes(x = x, y = y, col = growth), size = 0.1) + facet_wrap(~file)

        }

       
        if(dimension == 3){

            slice_thickness <- 10
            zslices <- unique((endStates %>% filter(trial == base))$zslice)
            sampleable_cells <- purrr::map(zslices, function(zslice){ inf_pos_pres %>% 
                filter(between(z, zslice - slice_thickness, zslice + slice_thickness)) %>% 
                 mutate(zslice = zslice) } ) %>% 
            bind_rows
            
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


            tps <- round(seq(2, mf, length.out = 5))
            toPlot <- left_join(sampleable_cells %>% filter(is.element(file, tps)),  inf_pos_pres_dth %>% select(-zslice) )

            edge <- toPlot %>%
                ggplot() + geom_point(aes(x = x, y = y, col = distToHull <= 20), size = 0.1) + facet_grid(zslice ~ file)
            pres <- toPlot %>%
                ggplot() + geom_point(aes(x = x, y = y, col = pressure < thres), size = 0.1) + facet_grid(zslice ~ file)
            drivs <- toPlot %>%
                ggplot() + geom_point(aes(x = x, y = y, col = growth), size = 0.1) + facet_grid(zslice ~ file)
           
        }

        relIndices <- (treeDat$maya_format)$index

        baseGrowth <- min(unique(inf_pos_pres_dth$growth))
        w <- as.numeric(gsub(".*_w|_d.*", "", base))

        suppInfForMayaFormat <- (inf_pos_pres_dth %>% ungroup() %>% filter(is.element(ID, relIndices)) %>% group_by(ID) %>% filter(file == min(file))  %>% 
            mutate(numDrivers = round(log(growth/baseGrowth)/log(w))) %>% select(ID, distFromCenter, distToHull, numDrivers) ) %>% ungroup() %>% rename(index = ID)

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

        write.table(maya_toPrint,
            file = paste0("../out/", outName, "/to_beast_format/", base, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

        avGrowth <- inf_pos_pres_dth %>% ungroup() %>% mutate(edgeAssociated = distToHull <= 20)  %>% 
            mutate(growthRate = growth * (pressure <= thres)) %>%
                group_by(file, edgeAssociated) %>% summarize(avGrowth = mean(growthRate)) %>%
                ggplot() + geom_line(aes(x = file, y = avGrowth, col = edgeAssociated)) +
                labs(x = "Time", y = "Average growth rate")

        write.table(
            inf_pos_pres_dth %>% ungroup() %>% mutate(edgeAssociated = distToHull <= 20)  %>% 
                mutate(growthRate = growth * (pressure <= thres)) %>%
                    group_by(file, edgeAssociated) %>% summarize(avGrowth = mean(growthRate), n = n()) %>% mutate(base = base),
            file = paste0("../out/", outName, "/growth_rate_differences/", base, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")


        figdir <- paste0("../out/", outName, "/figs/")

        system(paste0("mkdir ", figdir))
                
        pdf(paste0(figdir, "ancestral_locations_",base, ".pdf"), height = 5, width = 6.5 )
        print(plotGuesses(maya_toPrint))
        dev.off()
        
        pdf(paste0(figdir, "avGrowth_",base, ".pdf"), height = 5, width = 6 )
        print(avGrowth)
        dev.off()

        pdf(paste0(figdir, "inferredEdge_",base, ".pdf"), height = 8, width = 10 )
        print(edge)
        dev.off()

        pdf(paste0(figdir, "truePressure_",base, ".pdf"), height = 8, width = 10 )
        print(pres)
        dev.off()

        pdf(paste0(figdir, "drivers_",base, ".pdf"), height = 8, width = 10 )
        print(drivs)
        dev.off()


        if(dimension == 3){ treePlot <- plot_grid( g3, g4, g1, nrow = 1, rel_widths = c(1, 1, 1.5)) }
        if(dimension == 2){ treePlot <- plot_grid( plot_grid(g3, g4, nrow = 2), g1, nrow = 1) }
        
        pdf(paste0(figdir, "trees_",base, ".pdf"), height = 8, width = 8 )
        print(treePlot )
        dev.off()

    }
}


growth_rate_diffs <- foreach(fi = list.files(paste0("../out/", outName, "/growth_rate_differences/")), .combine = "rbind")%do%{

    tbl_df(read.table(
        file = paste0("../out/", outName, "/growth_rate_differences/", fi), header = TRUE, sep = "\t"))

}

growth_rate_diffs %>% separate(base, into = c("type", "m", "w", "d", "t", "mg", "mm", "l", "i", "s"), sep = "_") %>% 
    mutate(d = as.numeric(gsub("d", "", d))) %>% ggplot() + geom_point(aes(x = file, y = avGrowth  , col = edgeAssociated, group = edgeAssociated)) +
        geom_path(aes(x = file, y = avGrowth  , col = edgeAssociated, group = paste0(edgeAssociated, i))) +
            facet_wrap(~d, i) + labs(y = "Average birth rate", x = "Time")

pdf(paste0("../out/", outName, "/figs/growth_rate_differences_reps.pdf"), height = 5, width = 5)
print(growth_rate_diffs %>% separate(base, into = c("type", "m", "w", "d", "t", "mg", "mm", "l", "i", "s"), sep = "_") %>% 
    mutate(d = as.numeric(gsub("d", "", d))) %>% select(file, edgeAssociated, avGrowth, d, i, s) %>% spread(edgeAssociated, avGrowth) %>% 
        rename(center = `FALSE`, edge = `TRUE`) %>% mutate(diff = edge - center) %>% mutate(isna = is.na(diff), diff = ifelse(is.na(diff), 0, diff)) %>%
            filter(isna == FALSE, file < 110) %>% 
            ggplot() + geom_point(aes(x = file, y = diff/edge, col = factor(d)), alpha = 0.5) +
                geom_path(aes(x = file, y = diff/edge, col = factor(d), group = paste0(i, s)), alpha = 0.2) + labs(y =  "Birth rate difference (b_edge - b_center)/b_edge", x = "Time", col = "Death rate"))
dev.off()

pdf(paste0("../out/", outName, "/figs/growth_rate_differences.pdf"), height = 5, width = 5)
print(
    growth_rate_diffs %>% separate(base, into = c("type", "m", "w", "d", "t", "mg", "mm", "l", "i", "s"), sep = "_") %>% 
    mutate(d = as.numeric(gsub("d", "", d))) %>% select(file, edgeAssociated, avGrowth, d, i, s) %>% spread(edgeAssociated, avGrowth) %>% 
        rename(center = `FALSE`, edge = `TRUE`) %>% mutate(diff = edge - center) %>% mutate(isna = is.na(diff), diff = ifelse(is.na(diff), 0, diff)) %>%
            filter(isna == FALSE, file < 110)  %>% group_by(s, d) %>% filter(file > quantile(file, 0.1)) %>% summarize(diff = mean(diff/edge)) %>% ggplot() + geom_jitter(aes(x = d, y = diff, col = factor(d)), width = 0.01, height = 0) + labs(x = "Death rate", y = "Birth rate difference (b_edge - b_center)/b_edge", col = "death")
)
dev.off()

