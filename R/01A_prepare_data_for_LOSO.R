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