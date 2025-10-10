# curate full NCBI reference databases

# EvansAc_12S
source("R/settings_EvansAc_12S_ncbi_mammals.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_EvansAc_12S_ncbi_fishActinopterygii.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_EvansAc_12S_ncbi_fishNonActinopterygii.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_EvansAc_12S_ncbi_herps.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_EvansAc_12S_midori.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

# MiFishU_12S

source("R/settings_MiFishU_12S_ncbi_mammals.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_MiFishU_12S_ncbi_fishActinopterygii.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_MiFishU_12S_ncbi_fishNonActinopterygii.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_MiFishU_12S_ncbi_herps.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)

source("R/settings_MiFishU_12S_midori.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)


# Vences_16S

source("R/settings_Vences_16S_ncbi_mammals.R")
curate_amplicons(refdb = raw_refdb_path, fwd = fwd, rev = rev, taxon = NULL,
                 out = out_path, l = l, L = L, db_name = db_name, keep_all = F)
``