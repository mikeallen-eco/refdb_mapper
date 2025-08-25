# rules/15_blastg.smk

OUTPUT_DIR = config["paths"]["output"]

rule blastg:
    """
    Run blastg on each chunk of non-chimeric sequences.
    """
    input:
        chunk = f"{OUTPUT_DIR}/14_split/{{chunk}}.fasta"
    output:
        ghost_data = f"{OUTPUT_DIR}/15_blastg/{{chunk}}/ghost_data.csv",
        ghost_sum  = f"{OUTPUT_DIR}/15_blastg/{{chunk}}/ghost_sum.csv"
    conda:
        "envs/environment.yaml"
    params:
        reference_database = config["blastg"]["reference_database"],
        local_species = config["blastg"]["local_species"],
        marker = config["blastg"].get("marker", "average"),
        blast_path = config["blastg"].get("blast_path", "/usr/local/ncbi/blast/bin"),
        BLAST_args = config["blastg"].get("BLAST_args", "-max_target_seqs 5000 -perc_identity 85 -qcov_hsp_perc 90 -num_threads 4")
    run:
        import os

        # derive output folder from one of the output CSV paths
        outdir = os.path.dirname(output.ghost_data)
#        os.makedirs(params.outdir, exist_ok=True)

        shell(
            r"""
            mkdir -p "{OUTPUT_DIR}/15_blastg"
            
            blastg \
              --seqs="{input.chunk}" \
              --refdb="{params.reference_database}" \
              --locals="{params.local_species}" \
              --marker="{params.marker}" \
              --out="{outdir}" \
              --blast_path="{params.blast_path}" \
              --BLAST_args="{params.BLAST_args}" \
              --verbose
            """
        )
