# rules/14_split.smk

OUTPUT_DIR = config["paths"]["output"]
N_CHUNKS = config.get("n_tax_chunks", 10)
MIN_CLUSTER_SIZE = config.get("min_cluster_size", 2)
\
rule split:
    input:
        nochim_minsize = f"{OUTPUT_DIR}/13_min_size/pooled.cluster.consensus.nochim.min{MIN_CLUSTER_SIZE}.fasta"
    output:
        chunk = temp(f"{OUTPUT_DIR}/14_split/{{chunk}}.fasta")
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/14_split/{{chunk}}.log"
    params:
        outdir = f"{OUTPUT_DIR}/14_split"
    run:
        import os
        from math import ceil

        os.makedirs(params.outdir, exist_ok=True)

        # read sequences
        with open(input.nochim_minsize, "r") as f:
            lines = f.readlines()

        sequences = []
        header, seq_lines = None, []
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            sequences.append((header, "".join(seq_lines)))

        total_seqs = len(sequences)
        chunk_size = ceil(total_seqs / N_CHUNKS)

        # get chunk index from wildcard
        i = int(wildcards.chunk.split("_")[1]) - 1
        chunk_seqs = sequences[i*chunk_size:(i+1)*chunk_size]
        if not chunk_seqs:
            raise ValueError(f"No sequences for chunk {wildcards.chunk}")
        
        chunk_file = os.path.join(params.outdir, f"{wildcards.chunk}.fasta")
        with open(chunk_file, "w") as out_f:
            for h, s in chunk_seqs:
                out_f.write(f"{h}\n{s}\n")