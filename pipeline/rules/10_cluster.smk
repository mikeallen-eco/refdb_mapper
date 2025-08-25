# rules/10_cluster.smk

OUTPUT_DIR = config["paths"]["output"]

# clustering parameters
CLUSTER_IDENTITY = config.get("cluster_perc_id", 0.97)
VSEARCH_THREADS = config.get("vsearch_threads", 4)

rule cluster:
    input:
        pooled_fasta = f"{OUTPUT_DIR}/09_pool/pooled.fasta"
    output:
        clusters_centroids = temp(f"{OUTPUT_DIR}/10_cluster/pooled.cluster.centroids.fasta"),
        clusters_consensus = f"{OUTPUT_DIR}/10_cluster/pooled.cluster.consensus.fasta",
        clusters_uc = f"{OUTPUT_DIR}/10_cluster/pooled.cluster.uc",
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/10_cluster/10_cluster.log"
    shell:
        r"""
        mkdir -p {OUTPUT_DIR}/10_cluster

        vsearch \
            --cluster_fast {input.pooled_fasta} \
            --id {CLUSTER_IDENTITY} \
            --centroids {output.clusters_centroids} \
            --consout {output.clusters_consensus} \
            --uc {output.clusters_uc} \
            --sizein \
            --sizeout \
            --threads {VSEARCH_THREADS} \
            > {log} 2>&1

        echo "Clustering completed at identity {CLUSTER_IDENTITY}" >> {log}
        echo "Number of clusters: $(grep -c '^>' {output.clusters_centroids})" >> {log}
        """
