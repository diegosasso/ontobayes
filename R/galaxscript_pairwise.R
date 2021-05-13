#' Builds all input files and an script to run all pairwise analyses from ontobayes or ontobayes_profiles in Galax. 
#'
#' @param ids data.frame or list, if profile = FALSE, ids should be a data.frame providing matches between characters and anatomical entities. 
#'   The first column should provide character ids in the same order as presented in the character matrix. 
#'   The second column should provide the respective ontology ids referring to the anatomical entities. 
#'   The third column should provide the respective ontology terms referring to the anatomical entities.
#'   See examples provided in the data folder (data1.csv and data2.csv).
#'   Otherwise, ids should be a list of profiles organized as indicated in the examples.
#' @param profile logical, if TRUE specify that the input in ids is a list of profiles to query, organized as indicated in the examples.  
#' @param runs integer, the number of independent runs of mcmc used in the MrBayes commands block. Default is set to 2.
#' @param burnin integer, the absolute number of trees as the burn-in used in the mcmc. Default is set to 251.
#' @param OS character, sets the operational system to produce batch files for. Options are "windows" and "unix". Default is set to "windows".
#' @param filename character, name of the output file.
#'
#' @return All necessary input files, organized in subfolders inside the GALAX/INPUT/pairwise folder, 
#'   and an script to run all Galax analyses. 
#'
#' @examples
#' \dontrun{
#' # HAO example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("HAO:0000506",5), rep("HAO:0000513",5), rep("HAO:0000453",5), rep("HAO:0000234",5), 
#' rep("HAO:0001003",5), rep("HAO:0000874",5), rep("HAO:0000583",5), rep("HAO:0000630",5)),
#' c(rep("mandible",5), rep("maxilla",5), rep("labium",5), rep("cranium",5)
#' , rep("tentorium",5), rep("prothorax",5), rep("mesothorax",5), rep("metathorax",5))))
#'
#' # List of profiles with terms to query #
#' tags <- list(
#'   c("HAO:0000506", "HAO:0000513", "HAO:0000453"),
#'   c("HAO:0000234", "HAO:0001003"),
#'   c("HAO:0000874", "HAO:0000583", "HAO:0000630"))
#' 
#' # Set profile names and terms #
#' names(tags) <- c("mouthparts_profile1", "head_profile2", "mesosoma_profile3")
#' names(tags[[1]]) <- c("mandible", "maxilla", "labium")
#' names(tags[[2]]) <- c("cranium", "tentorium")
#' names(tags[[3]]) <- c("prothorax", "mesothorax", "metathorax")
#' 
#' # Creates bat file to run all analyses on Galax (profiles) #
#' galaxscript_pairwise(ids = tags, profile = T, runs = 2, burnin = 251, OS = "windows", filename = "galaxpairs1")
#'
#' # Creates bat file to run all analyses on Galax (all partitions) #
#' galaxscript_pairwise(ids = ID, profile = F, runs = 2, burnin = 251, OS = "windows", filename = "galaxpairsfull1")
#'
#' # UBERON example #
#' # Create ID object #
#' ID <- as.data.frame(cbind(paste0("C", 1:40),
#' c(rep("UBERON_0002244",5), rep("UBERON_0002397",5), rep("UBERON_0004742",5), rep("UBERON_2000658",5), 
#' rep("UBERON_2000488",5), rep("UBERON_0000151",5), rep("UBERON_0000152",5), rep("UBERON_0003097",5)),
#' c(rep("premaxilla",5), rep("maxilla",5), rep("dentary",5), rep("epibranchial bone",5)
#' , rep("ceratobranchial bone",5), rep("pectoral fin",5), rep("pelvic fin",5), rep("dorsal fin",5))))
#'
#' # List of profiles with terms to query #
#' tags <- list(
#'   c("UBERON_0002244", "UBERON_0002397", "UBERON_0004742"),
#'   c("UBERON_2000658", "UBERON_2000488"),
#'   c("UBERON_0000151", "UBERON_0000152", "UBERON_0003097"))
#' 
#' # Set profile names and terms #
#' names(tags) <- c("jaw_profile1", "pharynx_profile2", "fins_profile3")
#' names(tags[[1]]) <- c("premaxilla", "maxilla", "dentary")
#' names(tags[[2]]) <- c("epibranchial_bone", "ceratobranchial_bone")
#' names(tags[[3]]) <- c("pectoral_fin", "pelvic_fin", "dorsal_fin")
#' 
#' # Creates bat file to run all analyses on Galax (profiles) #
#' galaxscript_pairwise(ids = tags, profile = T, runs = 2, burnin = 251, OS = "windows", filename = "galaxpairs2")
#'
#' # Creates bat file to run all analyses on Galax (all partitions) #
#' galaxscript_pairwise(ids = ID, profile = F, runs = 2, burnin = 251, OS = "windows", filename = "galaxpairsfull2")
#' }
#' 
#' @export
galaxscript_pairwise <- function(ids = ids, profile = F, runs = 2, burnin = 251, OS = "windows", filename = "galaxscript")
{
  
  # Create input and output folders and subfolders for Galax #
  dir.create(file.path("GALAX", "INPUT",  "pairwise"), recursive = T)
  dir.create(file.path("GALAX", "OUTPUT", "pairwise"), recursive = T)
  
  if(profile == T){
    
    # Load all query terms #
    tags <- ids
    
    # Get all subfolder names #
    folders <- paste0(1:length(tags), "_", names(tags))
    
    x <- character()
    y <- character()
    
    for(i in 1:length(tags)){
      
      # Create a vector for the file names #
      x <- c(x, paste0(1:length(tags[[i]]), "_", names(tags[[i]])))  
      
      # Create a vector for the subfolder names #
      y  <- c(y, rep(folders[i], length(tags[[i]])))
      
    }
    
    if(OS == "windows"){
      
      # Create a vector with full path to the files #
      z <- paste0("NEXUS\\", y, "\\", x)
      
    }
    
    if(OS == "unix"){
      
      # Create a vector with full path to the files #
      z <- paste0("NEXUS/", y, "/", x)
      
    }
    
  }else{
    
    # Load all query terms #
    tags <- ids
    
    # Adjust ontology term names #
    tags <- gsub(tags[,3], pattern = " ", replacement = "_")
    tags <- gsub(tags, pattern = "/", replacement = "-")
	
    # Get all unique terms included in the dataset #
    x <- unique(tags)
    
    # Get all subfolder names #
    y <- paste0(1:length(x), "_", x)
    
    if(OS == "windows"){
      
      # Create a vector with full path to the files #
      z <- paste0("NEXUS\\ALL_PART\\", y, "\\", x)
      
    }
    
    if(OS == "unix"){
      
      # Create a vector with full path to the files #
      z <- paste0("NEXUS/ALL_PART/", y, "/", x)
      
    }	
    
  }
  
  # Create a vector for the names of all possible pairwise comparisons between terms #
  w <- combn(x, m = 2, simplify = F)
  
  # Create a vector for all possible pairwise comparisons between terms #
  v1 <- combn(z, m = 2, simplify = T)[1,]
  v2 <- combn(z, m = 2, simplify = T)[2,]
  
  block <- character()
  
  for(j in 1:length(v1)){
    
    for(k in 1:runs){
      
      # Write all pairwise comparison files #
      write(c(paste0(v1, ".nex.run", k, ".t")[j], paste0(v2, ".nex.run", k,".t")[j]), file = paste0("GALAX/INPUT/pairwise/", paste0(w[[j]], sep = "", collapse = "_"), ".run", k,".txt"))
      
      if(OS == "windows"){
        
        # Write all lines of pairwise comparisons in bat file #
        block <- c(block, paste0("galax -l GALAX\\INPUT\\pairwise\\", paste0(w[[j]], sep = "", collapse = "_"), ".run", k,".txt -s ", 
                                 burnin, " -g 1 -o GALAX\\OUTPUT\\pairwise\\galax_", paste0(w[[j]], sep = "", collapse = "_"), "_run", k)) 
        
      }
      
      if(OS == "unix"){
        
        # Write all lines of pairwise comparisons in sh file #
        block <- c(block, paste0("galax --listfile GALAX/INPUT/pairwise/", paste0(w[[j]], sep = "", collapse = "_"), ".run", k,".txt --skip ", 
                                 burnin, " --outgroup 1 --outfile GALAX/OUTPUT/pairwise/galax_", paste0(w[[j]], sep = "", collapse = "_"), "_run", k))
        
      }	
      
    }    
    
  }
  
  if(OS == "windows"){
    
    # Create headings for the bat file #
    headline <- paste(
      "@echo off",
      "set PATH=bin;%PATH%",
      "galax.exe %*",
      sep = "\n")
    
    # Write batch file for Windows #
    write(c(headline, block), file = paste0(filename, ".bat"))
    
  }
  
  if(OS == "unix"){
    
    # Create headings for the sh file #
    headline <- "#!/bin/bash"
    
    # Write batch file for Unix #
    write(c(headline, block), file = paste0(filename, ".sh"))
    
  }
  
}
