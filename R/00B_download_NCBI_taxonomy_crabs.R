# download the NCBI taxonomy files using crabs

download_NCBI_taxonomy_crabs <- function(out,
                                   conda_dir = "/Users/mikea/miniconda3/bin/conda",
                                   conda_env = "crb2"){
  
  # make working directory if needed
  if (!dir.exists(out)) {
    dir.create(out)
  }
  
  # make tax directory within out directory
  if (!dir.exists(paste0(out, "tax/"))) {
    dir.create(paste0(out, "tax/"))
  }
  
# download most recent NCBI tax files if needed
system2(conda_dir, 
        args = c("run", "-n", conda_env, "crabs", 
                 "--download-taxonomy",
                 "--output", paste0(out, "tax/")), 
        stdout = TRUE, stderr = TRUE)

}