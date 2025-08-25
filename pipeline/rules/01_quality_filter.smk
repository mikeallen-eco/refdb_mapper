RAW_DATA = config["paths"]["raw_data"]
OUTPUT_DIR = config["paths"]["output"]
QUALITY_SCORE = config.get("quality_score", 10)

rule quality_filter:
    input:
        f"{OUTPUT_DIR}/00_separate_samples/{{sample}}.fastq"
    output:
        qual = temp(f"{OUTPUT_DIR}/01_quality_filter/{{sample}}.q{QUALITY_SCORE}.fastq"),
        qual_distrib = f"{OUTPUT_DIR}/01_quality_filter/{{sample}}.q{QUALITY_SCORE}.distrib.txt"
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/01_quality_filter/logs/{{sample}}_quality_filer.log"

    shell:
        r"""
        echo "Quality filtering reads in {input} to >=q{QUALITY_SCORE}" >> {log}

        mkdir -p "{OUTPUT_DIR}/01_quality_filter"

        NanoFilt --q {QUALITY_SCORE} < {input} > {output.qual}

        chmod +x scripts/get_read_len_distribs.sh
        scripts/get_read_len_distribs.sh {output.qual} {output.qual_distrib}

        """
