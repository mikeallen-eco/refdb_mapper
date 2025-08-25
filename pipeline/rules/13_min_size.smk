# rules/13_min_size.smk

OUTPUT_DIR = config["paths"]["output"]
VSEARCH_THREADS = config.get("vsearch_threads", 4)
MIN_CLUSTER_SIZE = config.get("min_cluster_size", 2)

rule min_size:
    """
    Use seqkit to remove clusters with less than MIN_CLUSTER_SIZE reads.
    """
    input:
        nochim_fasta = f"{OUTPUT_DIR}/12_chimera/pooled.cluster.consensus.nochim.fasta"
    output:
        nochim_minsize = temp(f"{OUTPUT_DIR}/13_min_size/pooled.cluster.consensus.nochim.min{MIN_CLUSTER_SIZE}.fasta")
    conda:
        "envs/environment.yaml"
    params:
        min_cluster_size = config.get("min_cluster_size", 2)
    log:
        f"{OUTPUT_DIR}/13_min_size/13_min_size.log"
    shell:
        r"""
        mkdir -p {OUTPUT_DIR}/13_min_size
        
        ## Remove clusters under `min_cluster_size` reads

        ## generate list of numbers from 0 to min_cluster_size
        cluster_size_pattern=$(seq 0 $(({params.min_cluster_size} - 1)) | paste -sd'|' -)
        seqkit grep -r -n -v -p "size=($cluster_size_pattern)$" \
        "{input.nochim_fasta}" \
        -o "{output.nochim_minsize}"


        echo "Filtering based on minimum cluster size completed." >> {log}
        echo "Clusters smaller than {MIN_CLUSTER_SIZE} removed." >> {log}
        echo "$(grep -c '^>' {input.nochim_fasta}) clusters before removal."  >> {log}
        echo "$(grep -c '^>' {output.nochim_minsize}) clusters after removal."  >> {log}
        """