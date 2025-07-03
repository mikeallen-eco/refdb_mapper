# function to run leave-one-sequence-out loop for ghostblaster

# load libraries
library(dplyr)
library(tidyr)
library(Biostrings)
source("R/functions/LOOseq_refdb_subset.R")
source("R/functions/ghostblaster.R")
source("R/functions/subset_refdb.R")

LOOseq_ghostblaster <- function(refdb, # a reference database including the locals and more (e.g., all in genera or all regional species)
                                locals, # set of 'local' species to include in testing
                                ectoPredict,
                                out,
                                start_seq = 1,
                                BLAST_args = "-max_target_seqs 5000 -perc_identity 85 -qcov_hsp_perc 70",
                                verbose = F){
  # set arguments for testing
  testmode <- F
  if(testmode %in% T){
    refdb <- "data/midori_GB263_v16S.pga70.uni.l100.L300.n10.rdp.ama.fasta"
    locals <- read.csv("data/sp_names_tax.CdT20240712.csv")$orig_name
    # test_broader_db_only <- T
    out <- "data/LOOseq_GB_V16S_data/"
    ectoPredict = "data/ectoData_V16S/ectoTable/ectoTable_quantile.csv"
    BLAST_args = "-max_target_seqs 5000 -perc_identity 85 -qcov_hsp_perc 70"
    start_seq = 1
    verbose = T
    library(dplyr)
    library(tidyr)
    library(Biostrings)
    source("R/functions/LOOseq_refdb_subset.R")
    source("R/functions/subset_refdb.R")
  }
  
  # read in arguments to pass through to helper functions
  # needed?
  local_MOL_names <- locals
  GB_out <- out
  
  # create directories if needed
  if (!dir.exists(GB_out)) {
    dir.create(GB_out)
  }
  
  # make a temporary directory for diminished reference databases
  if (!dir.exists(paste0(GB_out, "tmp/"))) {
    dir.create(paste0(GB_out, "tmp/"))
  }
  
  # read in broader reference database
  r <- readDNAStringSet(refdb)
  
  # get df of broader reference database
  rn <- RDP_to_dataframe(r)
  
  # get refdb subset to Tumbira species 
  r.loc <- subset_refdb(RDPstringset = r,
                        locals = local_MOL_names)
  
  # get local species names that are in reference database (and have n>1)
  rn.loc <- RDP_to_dataframe(r.loc) %>%
    group_by(s) %>%
    mutate(n = length(s)) %>%
    ungroup() %>%
    filter(n > 1)
  sp_names.loc <- unique(rn.loc$s)
  
  # get vector of local sequences (with n > 1) to loop through in broader database
  local_seqnums_vector <- grep(paste(sp_names.loc, 
                                     collapse = "|"), names(r))
  
  # loop through species list to test for
  # allsp_results_list <- list()
  for (i in start_seq:length(local_seqnums_vector)) { 
    message("Testing sequence ", i, " of ", length(local_seqnums_vector), ": ", rn[local_seqnums_vector[i],]$s)
    tictoc::tic()
    
    seq_num <- local_seqnums_vector[i]
    
    # create diminished reference database subset (minus target sequence)
    refdb_dim <- LOOseq_refdb_subset(refdb, seq_num, return_db = TRUE)
    
    # write diminished reference database as a fasta in a tmp folder
    writeXStringSet(refdb_dim, 
                    filepath = paste0(GB_out,"tmp_GB_refdb.fasta"),
                    append = F)
    
    # create a reference database of just the LOO target sequence
    refdb_looseq <- LOOseq_refdb_subset(refdb, seq_num, return_db = FALSE) # , verbose = T
    
    # format the target sequence for taxonomy assignment
    seqs <- as.character(refdb_looseq)
    
    # run ghostblaster
    ghost_data <- ghostblaster(seqs,
                               refdb = paste0(GB_out, "tmp_GB_refdb.fasta"),
                               out = GB_out,
                               locals = locals,
                               ectoPredict = ectoPredict,
                               verbose = verbose,
                               BLAST_args = BLAST_args) 
    
    write.csv(ghost_data, 
              file = paste0(GB_out, "ghost_data_", 
                            ghost_data$qseqid[1],".csv"),
              row.names = F)
    
    tictoc::toc()
    
  }
  
}