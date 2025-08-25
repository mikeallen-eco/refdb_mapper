# rules/12_chimera.smk

OUTPUT_DIR = config["paths"]["output"]
VSEARCH_THREADS = config.get("vsearch_threads", 4)

rule chimera:
    """
    Detect and remove chimeric sequences from clustered consensus using vsearch (denovo method).
    """
    input:
        clusters_fasta = f"{OUTPUT_DIR}/10_cluster/pooled.cluster.consensus.fasta"
    output:
        nochim_fasta = temp(f"{OUTPUT_DIR}/12_chimera/pooled.cluster.consensus.nochim.fasta"),
        chimera_fasta = temp(f"{OUTPUT_DIR}/12_chimera/pooled.cluster.consensus.chimera.fasta"),
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/12_chimera/12_chimera.log",
    shell:
        r"""
        mkdir -p {OUTPUT_DIR}/12_chimera

        vsearch \
            --uchime_denovo {input.clusters_fasta} \
            --nonchimeras {output.nochim_fasta} \
            --chimeras {output.chimera_fasta} \
            --sizein \
            --sizeout \
            --threads {VSEARCH_THREADS} \
            > {log} 2>&1

        echo "Chimera check completed." >> {log}
        echo "Non-chimeric clusters (consensus): $(grep -c '^>' {output.nochim_fasta})" >> {log}
        echo "Chimeric clusters (consensus): $(grep -c '^>' {output.chimera_fasta})" >> {log}
        """
