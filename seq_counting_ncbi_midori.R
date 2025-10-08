
both <- RDP_to_dataframe(readDNAStringSet("~/Documents/mikedata/refdb_mapper/mammals_EvansAc_12S/refdb_EvansAc_12S_mammalia_ncbi20251008_midori265_tax20250609.fasta"), include_seqs = T)

ncbi <- RDP_to_dataframe(readDNAStringSet("~/Documents/mikedata/refdb_mapper/mammals_EvansAc_12S_ncbi/refdb_EvansAc_12S_mammalia_ncbi20251008_tax20250609.fasta"), include_seqs = T)

mid <- RDP_to_dataframe(readDNAStringSet("~/Documents/mikedata/refdb_mapper/mammals_EvansAc_12S_midori/refdb_EvansAc_12S_mammalia_midori265_tax20250609.fasta"), include_seqs = T)

length(unique(both$s))
length(unique(ncbi$s))
length(unique(mid$s))

length(unique(both$seqs))
length(unique(ncbi$seqs))
length(unique(mid$seqs))

length(unique(both$id))
length(unique(ncbi$id))
length(unique(mid$id))
