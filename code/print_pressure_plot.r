require(tidyverse)
require(foreach)
require(scatterplot3d)

dirpath <-  "../out/"

tumor_cols_base <- c("#F5E438", "#1D506E")
tumor_cols <- c("#F5E438", "#1D506E")
names(tumor_cols) <- c( TRUE, FALSE)
names(tumor_cols_base) <- c(TRUE, FALSE)


foreach(analysistype = c("2D_neut_bdg", "2D_sel_bdg_highermu", "3D_neut", "2D_sel_bdg", "3D_sel" ), .combine = "rbind")%do%{


    fullhullpath <- paste0(dirpath, analysistype,"/", "diversified_100/", "with_hull_info/")
    printpath <- paste0(dirpath, analysistype,"/", "pressure_plots/")
    system(paste0("mkdir ", printpath))

#    inds <- sort(c(seq(1, 225, by = 25)))
#can uncomment the above to just get a single example per d
    inds <- 1:225
    
    foreach(fi = list.files(fullhullpath)[inds], .combine = "rbind")%do%{

        hullinf <- tbl_df(read.table(paste0(fullhullpath, fi), header = TRUE))

        datToPlot <- hullinf %>% filter(file == max(file)) %>% arrange(desc(y))
        blueramp <- colorRampPalette(c(tumor_cols_base[1], "black"))
        maroonramp <- colorRampPalette(c(tumor_cols_base[2], "black"))       

        if(grepl("3D_neut", analysistype)){

            datToPlot <- datToPlot  %>% filter(x > 0)
            cols <- as.vector(tumor_cols[as.matrix(datToPlot %>% mutate(p = pressure > 1) %>% select(p)) + 1])
            rescale <- 1
            xlims <- range((datToPlot %>% select(x, y, z))$x)*rescale
            zlims <- range((datToPlot %>% select(x, y, z))$z)*rescale
            ylims <- range((datToPlot %>% select(x, y, z))$y)*rescale

            blues <- blueramp(ceiling(max((xlims)) * 1.5))

            maroons <- maroonramp(ceiling(max((xlims)) * 1.5))
            cexsize <- 0.5

            png(paste0(printpath, gsub("tsv", "png", fi)), height = 1.5, width = 1.5, units = "in", res = 300)
            par(mar = rep(0, 4))
            sc <- scatterplot3d(as.matrix(datToPlot %>% select(x, y, z)), pch = 20, color = cols, mar = c(0, 0, 0,0 ), cex.symbols =cexsize, angle = 0,
                                xlim = xlims, zlim = zlims, ylim = ylims, grid=FALSE, box=FALSE, axis = FALSE)
            stepsize <- -1
            for(v in seq(ceiling(max(xlims)), 0, by = stepsize)){
                dtp <- datToPlot %>% filter(between(x, v + stepsize, v))
                if(v <= 15){
                    tumor_cols_tmp <- c(blues[v + 1], maroons[v + 1])
                }else{
                    tumor_cols_tmp <- c(blues[v + 15], maroons[v + 15])
                }
                names(tumor_cols_tmp) <- names(tumor_cols_base)

                cols <- as.vector(tumor_cols_tmp[as.matrix(dtp %>% mutate(p = pressure > 1) %>% select(p)) + 1 ])
                sc$points3d(as.matrix(dtp %>% select(x, y, z)) , pch=16, col = cols, cex = cexsize)
            }
            dev.off()

        }

        if(grepl("3D_sel", analysistype)){

            datToPlot <- datToPlot  %>% filter(x > 0)
            cols <- as.vector(tumor_cols[as.matrix(datToPlot %>% mutate(p = pressure > 1) %>% select(p)) + 1])
            rescale <- 1
            xlims <- range((datToPlot %>% select(x, y, z))$x)*rescale
            zlims <- range((datToPlot %>% select(x, y, z))$z)*rescale
            ylims <- range((datToPlot %>% select(x, y, z))$y)*rescale

            blues <- blueramp(ceiling(max((xlims)) * 1.5))

            maroons <- maroonramp(ceiling(max((xlims)) * 1.5))
            cexsize <- 0.5

            png(paste0(printpath, gsub("tsv", "png", fi)), height = 1.5, width = 1.5, units = "in", res = 300)
            par(mar = rep(0, 4))
            sc <- scatterplot3d(as.matrix(datToPlot %>% select(x, y, z)), pch = 20, color = cols, mar = c(0, 0, 0,0 ), cex.symbols =cexsize, angle = 0,
                                xlim = xlims, zlim = zlims, ylim = ylims, grid=FALSE, box=FALSE, axis = FALSE)
            stepsize <- -1
            for(v in seq(ceiling(max(xlims)), 0, by = stepsize)){
                dtp <- datToPlot %>% filter(between(x, v + stepsize, v))
#dtp %>% ggplot() + geom_point(aes(x = x, y = y, col = numDriv)) 

                dtp <- dtp %>% mutate(numDriv = pmin(log(growth/min(hullinf$growth), base = 1.1), 4)) %>%
                        mutate(nd = paste(pressure <= 1, round(numDriv)))

                blues <- blueramp(5)
                maroons <- maroonramp(5)

                tumor_cols <- c(blues[1:5], maroons[1:5])
                names(tumor_cols) <- c(paste("TRUE", c(0:4)), paste("FALSE", c(0:4)))

                colkey <- bind_cols(tumor_cols, names(tumor_cols))
                names(colkey) <- c("col", "nd")
                
                cols <- as.vector((left_join(dtp, colkey) %>% select(col))$col)
                
                sc$points3d(as.matrix(dtp %>% select(x, y, z)) , pch=16, col = cols, cex = cexsize)
            }
            dev.off()

        }

        
        if(grepl("2D_neut", analysistype)){
            
            png(paste0(printpath, gsub("tsv", "png", fi)), height = 1.5, width = 1.5, units = "in", res = 300)
            print(datToPlot %>% arrange(ID) %>% ggplot() + geom_point(aes(x = x, y =y, col = pressure <= 1), size = 0.01) +
                scale_color_manual(values = tumor_cols_base) + theme_void() + theme(legend.position = "none")) + scale_size_continuous(range = c(0, 0.05))
            dev.off()
            
        }
        
        if(grepl("2D_sel_bdg", analysistype)){

            blues <- blueramp(5)
            maroons <- maroonramp(5)

            tumor_cols <- c(blues[1:5], maroons[1:5])
            names(tumor_cols) <- c(paste("TRUE", c(0:4)), paste("FALSE", c(0:4)))


#            names(tumor_cols) <- c("TRUE 0", "TRUE 1", "TRUE 2", "FALSE 0", "FALSE 1", "FALSE 2")
            toPlot <- (datToPlot %>% mutate(numDriv = pmin(log(growth/min(hullinf$growth), base = 1.1), 4)) %>% arrange(ID) %>%
                           mutate(nd = paste(pressure <= 1, round(numDriv))))


            png(paste0(printpath, gsub("tsv", "png", fi)), height = 1.5, width = 1.5, units = "in", res = 300)

            print(toPlot %>% 
                          ggplot() + geom_point(aes(x = x, y =y, col = nd), size = 0.01) +
                      scale_color_manual(values = tumor_cols) +
                          theme_void() + theme(legend.position = "none"))

            dev.off()
            
        }
    }
}




