#Ok, we want to regenerate the growth rate differences 


allRuns <- list.files("../out/")

foreach(i = allRuns)%do%{

    baseDir <- paste0("../out/", i, "/diversified_100/")
    system(paste0("mkdir ", baseDir, "updated_growth_rate_differences/"))

    tumorfiles <- list.files(paste0(baseDir, "with_hull_info/"))
    foreach(tum = tumorfiles)%do%{

        hullinf <- tbl_df(read.table(paste0(baseDir, "with_hull_info/", tum), header = TRUE))

        fi <- paste0(gsub("\\.tsv", "", tum), "_diversified_m1_n100.tsv")
        if(grepl("sigmoid", i)){
           
         write.table( hullinf %>% ungroup() %>% 
                mutate(edgeAssociated = distToHull <= 10)  %>% 
                mutate(growthRate = growth * (1 - 1/(1 +  exp(-5* (pressure - 1))))) %>%
                group_by(file, edgeAssociated) %>% 
                summarize(avGrowth = mean(growthRate), n = n()) %>% 
                mutate(base = base), file = fi,
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

        }else{
      
          write.table( hullinf %>% ungroup() %>% 
                mutate(edgeAssociated = distToHull <= 10)  %>% 
                mutate(growthRate = growth * (pressure <= thres)) %>%
                group_by(file, edgeAssociated) %>% 
                summarize(avGrowth = mean(growthRate), n = n()) %>% 
                mutate(base = base), file = fi,
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

        }
        
    }
}



