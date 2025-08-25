# rules/09_pool.smk

OUTPUT_DIR = config["paths"]["output"]

rule pool:
    """
    Pool per-sample dereplicated FASTAs into one file while retaining original sequence IDs.
    """
    input:
        lambda wildcards: [f"{OUTPUT_DIR}/08_dereplicate/{s}_dereplicate.fasta" for s in SAMPLES]
    output:
        temp(f"{OUTPUT_DIR}/09_pool/pooled.fasta")
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/09_pool/09_pool.log"
    shell:
        r"""
        mkdir -p {OUTPUT_DIR}/09_pool

        # Concatenate all per-sample dereplicated FASTAs
        cat {input} > {output}

        readcount=$(grep -o 'size=[0-9]*' {output} \
            | cut -d= -f2 \
            | paste -sd+ - \
            | bc)
        echo "Pooled FASTA contains ${{readcount}} total reads (accounting for ;size=X)." >> {log}

        echo "Pooled FASTA contains $(grep -c '^>' {output}) unique sequences." >> {log}
        """

