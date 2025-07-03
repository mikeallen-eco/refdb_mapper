# function that takes a MIDORI2 RDP-format reference database and extracts amplicons

extract_amplicons_crabs <- function(ref_seqs, fwd, rev, l, L,
                                   out){
  
  # make working directory if needed
  if (!dir.exists(out)) {
    dir.create(out)
  }
  
  if (verbose == TRUE) {message("Loading reference sequences into CRABS...")}
  # detect if refseqs is a fasta file and if not treat it like a DNAStringset
  if(grepl(ref_seqs[1], pattern = "fasta")){
    ref_seqs_path <- ref_seqs
  }else{
    writeXStringSet(ref_seqs, paste0(out, "raw_seqs.fasta"))
    ref_seqs_path <- paste0(out, "raw_seqs.fasta")
  }

if (verbose == TRUE) {message("Load reference sequences into CRABS and perform in silico PCR ...")}
tictoc::tic("Load reference sequences into CRABS and perform in silico PCR ...")
system2(conda_dir, 
        args = c("run", "-n", conda_env, "crabs", 
                 "--import",
                 "--import-format", "NCBI",
                 "--input", ref_seqs_path,
                 "--names", paste0(out, "tax/names.dmp"),
                 "--nodes", paste0(out, "tax/nodes.dmp"),
                 "--acc2tax", paste0(out, "tax/nucl_gb.accession2taxid"),
                 "--output", paste0(out, "crabs_raw_seqs.txt"),
                 "--ranks", "'superkingdom;phylum;class;order;family;genus;species'"
        ), 
        stdout = TRUE, stderr = TRUE)
tictoc::toc()

# In silico PCR to extract amplicons with primer regions (w/ ≤ 4 mismatches)
system2(conda_dir, 
        args = c("run", "-n", conda_env, "crabs", 
                 "--in-silico-pcr",
                 "--input", paste0(out, "crabs_raw_seqs.txt"),
                 "--output", paste0(out, "crabs_amplicons_pcr.txt"),
                 "--forward", fwd,
                 "--reverse", rev), 
        stdout = TRUE, stderr = TRUE)

tictoc::tic("Performing pga to extract amplicons without primer regions ...")
system2(conda_dir, 
        args = c("run", "-n", conda_env, "crabs", 
                 "--pairwise-global-alignment",
                 "--input", paste0(out, "crabs_raw_seqs.txt"),
                 "--amplicons", paste0(out, "crabs_amplicons_pcr.txt"),
                 "--output", paste0(out, "crabs_amplicons_pga.txt"),
                 "--forward", fwd,
                 "--reverse", rev,
                 "--size-select", "20000", 
                 "--percent-identity", "0.70",
                 "--coverage", "0.95"), 
        stdout = TRUE, stderr = TRUE)
tictoc::toc()

}
