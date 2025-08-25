#rules/length_filter.smk

import os

SAMPLES = config.get("samples", [])
RAW_DATA = os.path.expanduser(config["paths"]["raw_data"])
OUTPUT_DIR = os.path.expanduser(config["paths"]["output"])

QUALITY_SCORE = config.get("quality_score", 10)
MIN_LENGTH = config.get("raw_min_length", 100)
MAX_LENGTH = config.get("raw_max_length", 500)

rule length_filter:
    input:
        qfiltered_fastq = f"{OUTPUT_DIR}/01_quality_filter/{{sample}}.q{QUALITY_SCORE}.fastq"
    output:
        reads_len_filter = temp(f"{OUTPUT_DIR}/02_length_filter/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.fastq"),
        reads_len_filter_distrib = f"{OUTPUT_DIR}/02_length_filter/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.distrib.txt"    
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/02_length_filter/logs/{{sample}}_length_filter.log"

    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}/02_length_filter"
        
        set -euo pipefail
        input_fastq_no_ext=$(basename {input.qfiltered_fastq} .fastq)
        output_dir=$(dirname {output.reads_len_filter})

        echo "Length filtering reads in {input.qfiltered_fastq} between {MIN_LENGTH} and {MAX_LENGTH}" >> {log}

        seqkit seq --min-len {MIN_LENGTH} --max-len {MAX_LENGTH} {input.qfiltered_fastq} -o {output.reads_len_filter} >> {log} 2>&1
        
        chmod +x scripts/get_read_len_distribs.sh
        scripts/get_read_len_distribs.sh {output.reads_len_filter} {output.reads_len_filter_distrib} >> {log} 2>&1
        """

