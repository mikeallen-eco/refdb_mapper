# rules/11_otu_table.smk
OUTPUT_DIR = config["paths"]["output"]

rule otu_table:
    input:
        uc = f"{OUTPUT_DIR}/10_cluster/pooled.cluster.uc",
        centroids = f"{OUTPUT_DIR}/10_cluster/pooled.cluster.consensus.fasta"
    output:
        long = f"{OUTPUT_DIR}/11_otu_table/pooled.cluster.sample_counts.long.tsv",
        wide = f"{OUTPUT_DIR}/11_otu_table/pooled.cluster.otu_table.tsv"
    conda:
        "envs/environment.yaml"
    log:
        logfile = f"{OUTPUT_DIR}/11_otu_table/11_otu_table.log"
    run:
        import re, collections, os

        os.makedirs(f"{OUTPUT_DIR}/11_otu_table", exist_ok=True)

        def get_sample(label: str) -> str:
            """Extract sample name from label."""
            m = re.search(r"(?:^|[; ])sample=([^; ]+)", label)
            if m:
                return m.group(1)
            return re.split(r"[ _;]", label, maxsplit=1)[0]

        def get_size(label: str) -> int:
            """Extract read count from 'size=' in label; default to 1."""
            m = re.search(r"size=(\d+)", label)
            return int(m.group(1)) if m else 1

        def clean_label(label: str) -> str:
            """Remove size=... and centroid=... from label."""
            label = re.sub(r";size=\d+", "", label)
            label = re.sub(r"^centroid=", "", label)
            return label

        counts = collections.defaultdict(int)
        samples = set()
        centroids = set()

        with open(input.uc) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                rectype = parts[0]
                if len(parts) < 10:
                    continue  # malformed line

                query_label = parts[8]   # Col 9: original sequence from pooled file
                target_label = parts[9]  # Col 10: centroid label

                if rectype == "H":
                    sample = get_sample(query_label)
                    size = get_size(query_label)
                    centroid_id = clean_label(target_label)
                    counts[(centroid_id, sample)] += size
                    samples.add(sample)
                    centroids.add(centroid_id)

                elif rectype == "S":
                    sample = get_sample(query_label)
                    size = get_size(query_label)

                    if target_label == "*":
                        centroid_id = clean_label(query_label)
                    else:
                        centroid_id = clean_label(target_label)

                    counts[(centroid_id, sample)] += size
                    samples.add(sample)
                    centroids.add(centroid_id)

        samples = sorted(samples)
        centroids = sorted(centroids)

        # Long format
        with open(output.long, "w") as out:
            out.write("centroid\tsample\treads\n")
            for (centroid, sample), n in sorted(counts.items()):
                out.write(f"{centroid}\t{sample}\t{n}\n")

        # Wide (OTU-like) format
        with open(output.wide, "w") as out:
            out.write("centroid\t" + "\t".join(samples) + "\n")
            for centroid in centroids:
                row = [str(counts.get((centroid, sample), 0)) for sample in samples]
                out.write(centroid + "\t" + "\t".join(row) + "\n")

        with open(log.logfile, "w") as lg:
            lg.write(f"Centroids: {len(centroids)}\n")
            lg.write(f"Samples: {len(samples)}\n")
