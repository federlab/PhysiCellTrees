require(rdist)

sample_physicell_diversified_3D <- function(all_ids, tri,  n_sampled_cells = 100) {
    
    #filter for single tumor to sample 
    trial_ids <- all_ids %>% 
        dplyr::filter(trial == tri)

    num_z_slices <- 5
    slice_thickness <- 10
    quantiles <- seq(0.1, 0.9, length.out = num_z_slices)
    zslices <- (trial_ids %>% summarize(zpos = quantile(z, probs = quantiles)))$zpos

    #we have n_sampled_cells to work with in terms of our sampling
    #Let's figure out the number of cells at each of these z slices and then
    #sample numbers of cells proportionally

    sampleable_cells <- purrr::map(zslices, function(zslice){ trial_ids %>% 
           filter(between(z, zslice - slice_thickness, zslice + slice_thickness)) %>% 
           mutate(zslice = zslice) } ) %>% 
       bind_rows


    nslice <- sampleable_cells %>% group_by(zslice) %>% summarize(n = n()) %>% ungroup() %>%
        mutate(n_slice = rmultinom(1, n_sampled_cells, prob = n/n()))

    sample_helper <- function(trial_ids, n_sampled_cells){

    #randomize order of cell ids
        trial_ids <- trial_ids[sample(1:nrow(trial_ids), nrow(trial_ids), 
                                      replace = FALSE), ]
    
    print("pre pdist?")

    print("x:")
    print(trial_ids$x)
    print("y:")
    print(trial_ids$y)
    print("z:")
    print(trial_ids$z)

    print("lengths?")
    print(c(length(trial_ids$x), length(trial_ids$y), length(trial_ids$z)))

    tmp_mat <- as.matrix(data.frame(x = trial_ids$x, 
               			    y = trial_ids$y,
                                    z = trial_ids$z), ncol = 3)

    print(tail(tmp_mat))

    dist_mat <- rdist::pdist(tmp_mat)

    print("completed dist_mat step")
    #pairwise distance matrix between all alive cells in tumor
        dist_mat <- rdist::pdist(as.matrix(data.frame(x = trial_ids$x, 
                                                  y = trial_ids$y,
                                                  z = trial_ids$z), ncol = 3))

    print("pre fartherest point sampling?")

    #find sampled set of cells that maximized distance between sampled cells
        fps <- rdist::farthest_point_sampling(dist_mat, k = n_sampled_cells)
	

	print("post rdist?")
    
    #label cells as sampled or unsampled 
        trial_ids$sampled <- FALSE
        trial_ids$sampled[fps] <- TRUE
        
    #filter to only sampled cells 
        sampled_trial_ids <- trial_ids %>% 
            dplyr::filter(sampled)
        
        return(sampled_trial_ids)

    }

    purrr::map(zslices, function(zslice_in){
        numSamps <- (nslice %>% filter(zslice == zslice_in))$n_slice 
        sample_helper(sampleable_cells %>% filter(zslice == zslice_in), numSamps)
    }) %>% bind_rows
    
    
}



sample_physicell_diversified <- function(all_ids, tri, n_sampled_cells = 100) {
    
    #filter for single tumor to sample 
    trial_ids <- all_ids %>% 
        dplyr::filter(trial == tri)
    
    #randomize order of cell ids
    trial_ids <- trial_ids[sample(1:nrow(trial_ids), nrow(trial_ids), 
                                      replace = FALSE), ]
    
    #pairwise distance matrix between all alive cells in tumor
    dist_mat <- rdist::pdist(as.matrix(data.frame(x = trial_ids$x, 
                                                  y = trial_ids$y), ncol = 2))
    
    #find sampled set of cells that maximized distance between sampled cells
    fps <- rdist::farthest_point_sampling(dist_mat, k = n_sampled_cells)
    
    #label cells as sampled or unsampled 
    trial_ids$sampled <- FALSE
    trial_ids$sampled[fps] <- TRUE
    
    #filter to only sampled cells 
    sampled_trial_ids <- trial_ids %>% 
        dplyr::filter(sampled)
    
    return(sampled_trial_ids)
    
    
}

