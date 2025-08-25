RAW_DATA = config["paths"]["raw_data"]
OUTPUT_DIR = config["paths"]["output"]

rule separate_samples:
    input:
        pooled = RAW_DATA
    output:
        separated = temp(f"{OUTPUT_DIR}/00_separate_samples/{{sample}}.fastq"),
        separated_distrib = f"{OUTPUT_DIR}/00_separate_samples/{{sample}}.distrib.txt"
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/00_separate_samples/logs/{{sample}}_separate_samples.log"
    run:
        import re, gzip
        barcode_pattern = re.compile(r"barcode=(barcode\d+)")
        opener = gzip.open if input.pooled.endswith(".gz") else open

        with opener(input.pooled, "rt") as infile, open(output.separated, "w") as outfile:
            write_flag = False
            for line in infile:
                if line.startswith("@"):
                    m = barcode_pattern.search(line)
                    write_flag = (m and m.group(1) == wildcards.sample)
                if write_flag:
                    outfile.write(line)

        # length distrib script
        shell("""
            chmod +x scripts/get_read_len_distribs.sh
            scripts/get_read_len_distribs.sh {output.separated} {output.separated_distrib}
        """)
