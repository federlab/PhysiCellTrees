
setup_dirs <- function(dirname){
    
    system(paste0("mkdir ",dirpath,dirname))
    system(paste0("mkdir ",dirpath,dirname, "/config_files"))
    system(paste0("rm ",dirpath,dirname, "/config_files/*"))
    system(paste0("mkdir ",dirpath,dirname, "/raw"))
    system(paste0("mkdir ",dirpath,dirname, "/processed"))
    system(paste0("mkdir ",dirpath,dirname, "/fig"))

}

runner_restarter <- function(dirname, command = 'tumor_evolve'){

    prefix <- paste0(dirpath,dirname,'/config_files/')
    confsToRun <- (list.files(paste0(prefix)))

    registerDoParallel(cores = detectCores() )

    foreach(filename = confsToRun)%dopar%{

        configpath <- paste0(prefix, filename)
        print(paste0("./", command, " ",configpath, " "))
        system(paste0("/net/feder/vol1/project/tumor_evolutionary_sims/PhysiCell/", command, " ",configpath, " "))

        #We here want to check to see that the job finished with a non-0 number of cells
        #We'll do this by checking the file size

	endsOn0State = TRUE

	finalState_file <- paste0(dirpath,dirname, "/raw/", gsub("\\.xml", "", filename), "/final_cells_physicell.mat")
	if(file.exists(finalState_file)  & file.info(finalState_file)$size > 100){
                endsOn0State = FALSE
	}

        while(endsOn0State){
            #restart the job with a different seed
            
            #So, the random seed is written into the xml. We'll need to edit. 
            xmlin <- read_file(file = paste0(prefix, filename))
            xmlout <- gsub( "[0-9]*</random_seed>", paste0(sample(1:100000, 1), "</random_seed>"), xmlin)
            write_file(x = xmlout, file = paste0(prefix, filename))

            #Let's clear the raw output folder
            system((paste0("rm ", dirpath,dirname, "/raw/", gsub("\\.xml", "", filename), "/*")))
            
            #Then, we restart the task
            system(paste0("/net/feder/vol1/project/tumor_evolutionary_sims/PhysiCell/", command, " ",configpath, " "))

            #When the task ends, we see if endsOn0State has changed
             if(file.exists(finalState_file) & file.info(paste0(dirpath,dirname, "/raw/", gsub("\\.xml", "", filename), "/final_cells_physicell.mat"))$size > 100 ){
                 endsOn0State = FALSE
             }
        }

    }

}


process <- function(dirname){

    dirsh <- paste0(dirpath,dirname,"/")
    dpath <- paste0(dirsh, "raw/")

    foreach(dirv = list.files(dpath))%do%{

        filestoparse <- list.files(paste0(dpath, dirv))
        filestoparse <-  filestoparse[grep("output.*_physicell\\.mat", filestoparse)]

        maxi <- max(as.numeric(gsub("_.*|[A-Za-z]|\\.", "", 
                                    filestoparse)), na.rm = TRUE)

        if(file.exists(paste0(dpath, dirv, "/final_cells_physicell.mat"))){ filestoparse <- c(filestoparse, "final_cells_physicell.mat") }				    

        alldat <- foreach (i = filestoparse, .combine = "rbind")%do%{

	    indat <- (t(readMat(paste0(dpath, dirv, "/", i))$cells))

            colnames(indat) <- c("ID", "position.x", "position.y", "position.z",
                         "total_volume", "cell_type", "cycle_model", "current_phase",
                         "elapsed_time_in_phase", "nuclear_volume", "cytoplasmic_volume",
                         "fluid_fraction", "calcified_fraction", "orientation.x",
                         "orientation.y", "orientation.z", "polarity", "migration_speed",
                         "motility_vector.x", "motility_vector.y", "motility_vector.z",
                         "migration_bias", "motility_bias_direction.x",
                         "motility_bias_direction.y", "motility_bias_direction.z",
                         "persistence_time", "motility_reserved", "mutated", 
                         "growth", "mig_mut",  "mig", "pressure")
            indat <- tbl_df(indat)

            pr <- indat %>% 
                mutate(file = as.numeric(gsub("[A-Za-z]|\\.|_", "", i))) %>%
                mutate(file = ifelse(is.na(file), maxi + 1, file)) %>% 
                mutate(x = round(position.x, 3), 
                       y = round(position.y,3), 
                       z = round(position.z, 3)) %>%
		mutate(pressure = round(pressure, 3)) %>%
                mutate(rad = (3/(4*pi)*total_volume)^(1/3)) %>% group_by(ID) %>% 
                ungroup() %>%
                dplyr::select(ID, x, y, file, pressure, mutated, growth, mig_mut, mig, z)

        }

        write.csv(alldat , file =
			 paste0(dirsh, "processed/", dirv, ".csv"), quote = FALSE)

    }
}


parse_relationships <- function(dirname){

    dirsh <- paste0(dirpath,dirname,"/")
    dpath <- paste0(dirsh, "raw/")

    registerDoParallel(cores = detectCores() )

    foreach(dirv = list.files(dpath))%dopar%{

        logfile <- paste0(dpath, dirv, "/log.csv")
        relationshipTree <- rawToDescGraph(logfile)

        write.csv(relationshipTree$array_format , file =
                      paste0(dirsh, "processed/", dirv, "_par_child.csv"), quote = FALSE)

        names(relationshipTree$child_list) <- paste(1:length(relationshipTree$child_list))

        array_of_relationship_tree <- do.call("rbind", mapply(function(x, i){
            if(is.null(x)){ x = "NULL" }
            cbind(i, x)
          }, relationshipTree$child_list, names(relationshipTree$child_list)))


        write.csv(array_of_relationship_tree, file =
                  paste0(dirsh, "processed/", dirv, "_child_list.csv"))
        
    }
}

setup_configs_cube <- function(width = 1000, height = 1000, depth = 1000, threads = 4, migration = 0, 
                               fitness = 1, death = 0, threshold = 10, 
                          mutation_growth = 1, mutation_mig = 1, time = 5000, reps = 3,
			  		      cellNumEndCondition = 10000, driverLim = 3, initialCellNumber = 1,
                          dirname, use2D = 'true'){

    setup_dirs(dirname)

    enable_z = 0
    if(use2D != "true"){
        enable_z = 1
    }

    params <- expand.grid(migration, fitness, death,
                          threshold, mutation_growth, mutation_mig, 
                          time, cellNumEndCondition, driverLim, initialCellNumber)
    names(params) <- c("migration", "fitness", "death",
                       "threshold", "mutation_growth", "mutation_mig", 
                       "time", "cellNumEndCondition", "driverLim", "initial_cell_number")

    foreach(par = 1:nrow(params))%do%{
        foreach(repval = 1:reps)%do%{

    	    rseed <-  sample(1:100000, 1)

            p <- params[par,]
        
            pstring <- paste0("sampconfig_m",p$migration,
                              "_w",p$fitness,
                              "_d",p$death,
                              "_t",p$threshold,
                              "_mg",p$mutation_growth,
                              "_mm",p$mutation_mig,
                              "_l",p$time,
                              "_i",repval, 
						      "_s", rseed) 

            relDir <- paste0(dirpath,dirname, "/raw/", pstring)
            
            xmltoprint <- sprintf('<?xml version="1.0" encoding="UTF-8"?>

<PhysiCell_settings version="devel-version">
		    <domain>
			<x_min>-%d</x_min>
				<x_max>%d</x_max>
					<y_min>-%d</y_min>
						<y_max>%d</y_max>
							<z_min>-%d</z_min>
								<z_max>%d</z_max>
									<dx>20</dx>
										<dy>20</dy>
											<dz>20</dz>
												<use_2D>%s</use_2D>
												</domain>
												
												<overall>
													<max_time units="min">%d</max_time> 
														  <time_units>min</time_units>
															<space_units>micron</space_units>
															</overall>
															
															<parallel>
															 <omp_num_threads>%d</omp_num_threads>
															 </parallel> 
															 
															 <save>
															  <folder>%s</folder> <!-- use . for root --> 

															  		      <full_data>
																	        <interval units="min">360</interval>
																			   <enable>true</enable>
																			    </full_data>
																			     
																			      <SVG>
																			        <interval units="min">60</interval>
																					   <enable>false</enable>
</SVG>
		
			<legacy_data>
					<enable>false</enable>
						</legacy_data>
						</save>

<options>
		<virtual_wall_at_domain_edge>true</virtual_wall_at_domain_edge>
<legacy_random_points_on_sphere_in_divide>false</legacy_random_points_on_sphere_in_divide>

	</options>

	<user_parameters>
		<random_seed type="int" units="dimensionless">%d</random_seed> 
			     <cell_migration_speed type="double" units="micron/min">%f</cell_migration_speed>	
			     			   <cell_relative_cycle_entry_rate type="double" units="dimensionless">%f</cell_relative_cycle_entry_rate>
						   				   <cell_death_rate type="double" units="dimensionless">%f</cell_death_rate>
										   		    <pressure_threshold type="double" units="dimensionless">%f</pressure_threshold>
												    			<cell_mutation_rate_growth type="double" units="dimensionless">%f</cell_mutation_rate_growth>
																		   <cell_mutation_rate_mig type="double" units="dimensionless">%f</cell_mutation_rate_mig>
                <cellNumEndCondition type="double" units="dimensionless">%f</cellNumEndCondition>
                <initial_cell_number type="double" units="dimensionless">%f</initial_cell_number>
                <driverLim type="double" units="dimensionless">%f</driverLim>
                <print_directory type="string" units="dimensionless">%s</print_directory>
                <enable_z type="int" units="dimensionless">%s</enable_z>
		</user_parameters>
</PhysiCell_settings>)',  width, width, height, height, depth, depth, use2D, time, threads,
                          relDir, rseed,
                          p$migration, p$fitness, 
                          p$death, p$threshold, p$mutation_growth, p$mutation_mig, 
			  	   		  p$cellNumEndCondition, p$initial_cell_number, p$driverLim, 
                          relDir, enable_z)

  

            write(xmltoprint, paste0(dirpath, dirname, "/config_files/", pstring, ".xml"))
            system(paste0("mkdir ",dirpath,dirname, "/raw/", pstring))
        }
    }
}
   

