require(foreach)
require(tidyverse)
require(viridis)
require(rgl)
require(phytools)
require(ggtree)
require(cowplot)
require(R.matlab)


#This is a slow step of figuring out which cells are descendants of which
#I split the logic from these functions in case we want to resample

#Function: rawToDescGraph
#
#Inputs:
#Path to a log file 'directory' containing a tumor history
#
#Outputs:
#A list containing two entries:
# 1) array_format - spatial position data for all cells
# 2) child_list - relationship data for all cells

rawToDescGraph <- function(directory, numCellsToSample = 20, mu = 1, rseed = NA, samplingScheme = NA){


#Read in log files printed by PhysiCell
    log <- read.table(directory, sep = ",")
    colnames(log) <- c("t", "par", "x.par", "y.par", "z.par",
                       "daugh1", "x.1", "y.1", "z.1",
                       "daugh2", "x.2", "y.2", "z.2", "rad")
    log <- tbl_df(log)

#this seems to be the easiest way to gather here
#Each PC print statement lists both daughters - we want each daughter to have its own row
    d1 <- log %>% dplyr::select(t, par, x.par, y.par, z.par, daugh1, x.1, y.1, z.1) %>%
        rename(ID = daugh1, x = x.1, y = y.1, z = z.1)
    d2 <- log %>% dplyr::select(t, par, x.par, y.par, z.par, daugh2, x.2, y.2, z.2) %>%
        rename(ID = daugh2, x = x.2, y = y.2, z = z.2)
    dat <- bind_rows(d1, d2) %>% arrange(t, par, ID)

#Mirror Maya's file format type and names
    maya_format <- left_join(dat, 
          dat %>% dplyr::select(t, par) %>% unique() %>% rename(ID = par, deathdate = t)) %>% 
              rename(birthdate = t, index = ID, locx = x, locy = y, locz = z, 
              parent_index = par) %>% 
              dplyr::select(birthdate, deathdate, index, locx, locy, locz, parent_index) %>% 
              arrange(index)

#Ok, now we need a way to allow each birth event to permit a mutation with a certain
#probability

#For each index, let's make a list of descendants
#we can start at the bottom so we save some work

    child_list <- list()
    maxind <- max(maya_format$index)
    child_list[[maxind]] <- NA

    for(ind in sort(unique(maya_format$index), decreasing = TRUE)){
        if(ind %% 5000 == 0){print(paste0(ind, " nodes remaining"))}
#All nodes are their own descendants
        mychildren <- (maya_format %>% filter(parent_index == ind))$index
        child_list[[ind]] <- c(child_list[[mychildren[1]]], child_list[[mychildren[2]]], ind)
    }

    return(list(array_format = maya_format, 
                child_list = child_list))
}





#Function: sampleTree
#Inputs:
#Path to a log file containing a tumor history
#Sampling scheme - "random", "grid" (default: random)
#Sample size (default = 20)
#Mutation rate
#Random seed

#Outputs:
#A list containing four entries:
# 1) toBD - for piping into the Brownian Diffusion analysis (tbl)
# 2) timetree - the TRUE time tree (newick)
# 3) muttree - the true mutation tree (newick)
# 4) maya_format - this goes to Maya



sampleTree <- function(relationshipTree, numCellsToSample = 20, mu = 1, rseed = NA, samplingScheme = NA, sampleIDs){
rpois( length((relationshipTree$array_format)$index), 1)

    maya_format = relationshipTree$array_format
    child_list = relationshipTree$child_list
    
    if(is.na(rseed)){
        rseed <- floor(runif(1, 1, 1000000))
    }

#births_to_mark <- rep( (relationshipTree$array_format)$index, (rpois( length((relationshipTree$array_format)$index), mu)))
    births_to_mark <- rep( names(child_list), (rpois(length(child_list), mu)))

#Let's sort them so mutations occur in order
    births_to_mark <- sort(births_to_mark)
#Most of these will be late/low frequency
    
#Muts is going to be a text string that we construct to keep track of which mutations
# are present in any given cell

#What if we made an intermediate index that moves indices to row numbers?

    muts <- rep("", max(maya_format$index))
    for(mut_num in 1:length(births_to_mark)){
    	#Old (potentially problematic version)
	#        muts[child_list[[births_to_mark[mut_num]]] ] <- 
	#            paste0(muts[child_list[[births_to_mark[mut_num]]] ], mut_num, sep = ",")

        muts[child_list[[paste(births_to_mark[mut_num])]] ] <- 
            paste0(muts[child_list[[paste(births_to_mark[mut_num])]] ], mut_num, sep = ",")

    }
    muts <- muts[maya_format$index]


#Let's filter out the last comma (where present) and add brackets
# to match Maya's format
    muts_formatted <- paste0("[",(gsub(",$", "", muts)), "]")
    maya_format_withmuts <- maya_format %>% 
        mutate(mutations = muts_formatted) %>% 
            mutate(numMuts = str_count(mutations, "[0-9][,\\]]")) 

#Subtle but perhaps important point that's different between PC sims and Maya's sims - 
# Where these cells are born is NOT where they're sampled, necessarily. In this log file, 
# I'm only tracking where things are born. I need to actually merge in where things are sampled
# at some point. This info is also saved but in a different file that I need to dig up separately. 

#Ok, let's turn these sampled cells into a tree
#Should set a seed
    set.seed(rseed)
    toPlot  <-  maya_format_withmuts
    sampledLiveCells <- toPlot %>% #filter(is.na(deathdate)) %>% 
#        filter(between(locz, -20, 20)) %>% 
        filter(is.element(index, sampleIDs)) %>% 
        sample_n(numCellsToSample) 

#allAnc (i.e., all ancestors) will keep track of ALL divisions that lead to extant cells
#This could be helpful for reconstructing lineage trajectories, but has more info
# than what we want for our trees. We'll write something to trim extra cell divisions later
    allAnc <- sampledLiveCells %>% mutate(extant = TRUE)
    for(ind in sampledLiveCells$index){
        print(paste0(which(ind == sampledLiveCells$index), "/", 
                     length(sampledLiveCells$index)))
        currind <- ind
        parind <- (allAnc %>% filter(index == currind))$parent_index
    #As long as your parent isn't in the plot set, add it
    #This keeps us from adding everything a bunch of times
        while(!is.element(parind, allAnc$index) & parind > 0){
            allAnc <- allAnc %>% 
                add_row(toPlot %>% filter(index == parind) %>% mutate(extant = FALSE))
            currind <- parind
            parind <- (allAnc %>% filter(index == currind))$parent_index
        }
    }


#Also, I think we need to add a 0 entry because it's missing because it technically was never born

    allAnc <- allAnc %>% add_row(birthdate = 0, 
                    deathdate = 0, #this isn't correct, but ok for now
                    index = 0, 
                    locx = 0,  #I think I actually don't necessarily initiate at 0,0,0
                    locy = 0,  #(it's randomly dist around 0,0,0)
                    locz = 0,  # but I can change that in future sims
                    parent_index = NA, 
                    mutations = "[]", 
                    numMuts = 0, 
                    extant = FALSE) %>% 
                    arrange(birthdate)

#I think we would like a version where we cut out intermediate relationships for cells
#that are not MRCAs 
#Let's first identify which cells are the MRCAs (they should be the only ones with two children)
    mrcas <- (allAnc %>% group_by(parent_index) %>% 
                  summarize(n = n()) %>% filter(n == 2))$parent_index

#So, our goal is to create a version in which cells are only
#included if they're extant or MRCAs
    allTree <- bind_rows(allAnc %>% filter(extant == TRUE),
                         allAnc %>% filter(is.element(index, mrcas)))

    for(ind in allTree$index){

        print(paste0(which(ind == allTree$index), "/", 
                     length(allTree$index)))
        currind <- ind
        parind <- (allTree %>% filter(index == currind))$parent_index

    #If a cell's listed parent is not an mrca, keep moving through the file
        while(!is.element(parind, mrcas) & !is.na(parind)){
            currind <- parind
            parind <- (allAnc %>% filter(index == currind))$parent_index
        }

    #When we reach an MRCA, update the parent
        allTree <- allTree %>% 
            mutate(parent_index = ifelse(index == ind, parind, parent_index))
    }


#Ok, at this point, this should by what we need for tree plotting
#Start with the node that's most ancestral

    makeNewick_muttree <- function(allTree){

        currnode <- root

        recurse_helper <- function(currnode){

            kids <- allTree %>% filter(parent_index == currnode$index)

            if(nrow(kids) == 0){
                return(currnode$index)
            }

            kid1 <- recurse_helper(kids[1,])
            kid2 <- recurse_helper(kids[2,])

            kidmuts <- kids$numMuts - currnode$numMuts

            paste0("(", kid1, ":",kidmuts[1], ",", kid2, ":", kidmuts[2], ")", currnode$index)
        }

        phylab <- recurse_helper(root)
        read.newick(text = paste0(phylab, ";\n"))

    }

    makeNewick_timetree <- function(aT){

        currnode <- root

        recurse_helper <- function(currnode, aT){

            kids <- aT %>% filter(parent_index == currnode$index)

            if(nrow(kids) == 0){
                return(currnode$index)
            }

            kid1 <- recurse_helper(kids[1,], aT)
            kid2 <- recurse_helper(kids[2,], aT)

            kidtimes <- kids$deathdate - currnode$deathdate

            paste0("(", kid1, ":",kidtimes[1], ",", kid2, ":", kidtimes[2], ")", currnode$index)
        }

        phylab <- recurse_helper(root, aT)

        read.newick(text = paste0(phylab, ";\n"))

    }


    root <- allTree %>% filter(is.na(parent_index))

    nwk <- makeNewick_muttree(allTree)

    guessdeath <- max(allTree$birthdate) + 4
#How do we make a time tree?
    allTree_wdeath <- allTree %>% mutate(deathdate = ifelse(is.na(deathdate), guessdeath, deathdate))

    nwk_tt <- makeNewick_timetree(allTree_wdeath)


#What do we still need here? 
#branch lengths (in numbers of mutations)
#and time differences (in terms of minutes)
#we'll need to estimate the time differences based on branch lengths

    parMutsAndTime <- allTree_wdeath %>% dplyr::select(index, birthdate, numMuts) %>% 
        rename(parent_index = index, parbirth = birthdate, parmuts = numMuts)

    toBD <- left_join(allTree_wdeath, parMutsAndTime) %>% 
        mutate(b = numMuts - parmuts, deltat = birthdate - parbirth) %>% 
            dplyr::select(-c(parbirth, parmuts, mutations, deathdate)) %>% 
                rename(x = locx, y = locy, z = locz, id = index, par = parent_index, 
                       fixed = extant, kup = numMuts) %>% 
                           dplyr::select(id, par, x, y, z, fixed, b, deltat, kup, birthdate )

    #This script can sort of go into a few different formats. 
    #For now, I'll return each of them as an element of a list, but maybe in the future
    #I'll make bespoke printing


    return(list(toBD = toBD, 
                timetree = nwk_tt,
                muttree = nwk, 
                maya_format = allAnc))

}


