# rules/17_join_otu_tables.smk

OUTPUT_DIR = config["paths"]["output"]

rule join_otu_tables:
    """
    Join ghost_sum.csv (taxonomy assignment by OTU) with OTU read counts created in otu_reads_by_sample rule.
    """
    input:
        combined_sum  = f"{OUTPUT_DIR}/16_combine_blastg/ghost_sum.csv",
        otu_table = f"{OUTPUT_DIR}/11_otu_table/pooled.cluster.otu_table.tsv"
    output:
        joined_otus = f"{OUTPUT_DIR}/17_join_otu_tables/final_otu_table.csv",
    conda:
        "envs/environment.yaml"
    params:
        outdir = f"{OUTPUT_DIR}/17_join_otu_tables"
    log:
        joined_otus_log = f"{OUTPUT_DIR}/17_join_otu_tables/join_otu_tables.log"
    shell:
        r"""
        mkdir -p "{params.outdir}"
        
        Rscript scripts/join_otu_tables.R \
            --otu_table_tsv {input.otu_table} \
            --ghost_sum_csv {input.combined_sum} \
            --out_file "{output.joined_otus}"

        num_otus=$(( $(wc -l < {output.joined_otus}) - 1 ))
        echo "Join complete." >> {log.joined_otus_log}
        echo "Writing CSV with read count and taxonomic assignment info for ${{num_otus}} clusters" >> {log.joined_otus_log}

        """
