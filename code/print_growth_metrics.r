#Compute the growth rate of edge and center-associated cells in physicell tumors

computeGrowthRates <- function(dirToProcess, sampleSchemeName){


    outName <- sampleSchemeName
    reldir <- paste0("../out/", dirToProcess, sampleSchemeName, "/with_hull_info/")

    system(paste0("mkdir ../out/", dirToProcess, sampleSchemeName,"/growth_rate_differences/"))

    list_of_dirs <- list.files(reldir)

    foreach(dirstub = list_of_dirs)%do%{

        inf_pos_pres <- tbl_df(read.table(paste0(reldir, dirstub),header = TRUE, sep = "\t"))
        base <- gsub("\\.tsv", "", dirstub)
        
        thres <- as.numeric(gsub(".*_t|_mg.*", "", dirstub))

        write.table(
            inf_pos_pres %>% ungroup() %>% mutate(edgeAssociated = distToHull <= 10)  %>% 
                mutate(growthRate = growth * (pressure <= thres)) %>%
                group_by(file, edgeAssociated) %>% summarize(avGrowth = mean(growthRate), n = n()) %>% 
                    mutate(base = base ),
            file = paste0("../out/", dirToProcess, "/", outName, "/growth_rate_differences/", base, ".tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
                
    }
}


computeGrowthRates("2D_neut_bdg/", "diversified_100")
computeGrowthRates("2D_sel_bdg_highermu/", "diversified_100")
computeGrowthRates("3D_neut/", "diversified_100")


#Computer the total number of edge or center associated cells at the sampling stage of a physicell tumor

edgeVCenterCounts <- function(dirToProcess, sampleSchemeName){

    outName <- sampleSchemeName
    reldir <- paste0("../out/", dirToProcess, sampleSchemeName, "/with_hull_info/")

    list_of_dirs <- list.files(reldir)

    runscheme <- gsub(".*/", "", gsub("/$", "", dirToProcess))

    edgecounts <- foreach(dirstub = list_of_dirs, .combine = "rbind")%do%{

        inf_pos_pres <- tbl_df(read.table(paste0(reldir, dirstub),header = TRUE, sep = "\t"))

        inf_pos_pres %>% filter(file == max(file)) %>% 
            mutate(edge = ifelse(distToHull <= 20, 'edge', 'center')) %>% 
            group_by(edge) %>% count() %>% mutate(file = dirstub) %>% spread(edge, n)  %>% 
            mutate(runScheme = runscheme)
        
    }


    write.table(edgecounts, file = paste0("../out/", dirToProcess, "edge_v_center_counts_", runscheme, ".tsv"),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
}


edgeVCenterCounts("2D_neut_bdg/", "diversified_100")
edgeVCenterCounts("2D_sel_bdg_highermu/", "diversified_100")
edgeVCenterCounts("3D_neut/", "diversified_100")


