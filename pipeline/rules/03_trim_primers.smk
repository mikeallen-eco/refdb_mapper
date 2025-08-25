# rules/03_trim_primers.smk
QUALITY_SCORE = config.get("quality_score", 10)
MIN_LENGTH = config.get("raw_min_length", 100)
MAX_LENGTH = config.get("raw_max_length", 500)
OUTPUT_DIR = config["paths"]["output"]
CUTADAPT_ERROR = config.get("cutadapt_error_rate", 0.2)

rule trim_primers:
    input:
        reads_len_filter = f"{OUTPUT_DIR}/02_length_filter/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.fastq",
        primer_yaml = "config/pipeline_config.yaml"
    output:
        linked = temp(f"{OUTPUT_DIR}/03_trim_primers/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked.fastq"),
        linked_dd = temp(f"{OUTPUT_DIR}/03_trim_primers/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked.dd.fastq"),
        unlinked = temp(f"{OUTPUT_DIR}/03_trim_primers/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.unlinked.fastq"),
        unlinked_sub = temp(f"{OUTPUT_DIR}/03_trim_primers/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.unlinked.sub.fastq"),
        linked_dd_distrib = f"{OUTPUT_DIR}/03_trim_primers/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked.dd.distrib.txt"
    conda:
        "envs/environment.yaml"
    log:
        trim_primers_log = f"{OUTPUT_DIR}/03_trim_primers/logs/{{sample}}_trim_primers.log"
    shell:
        r"""
        mkdir -p "{OUTPUT_DIR}/03_trim_primers"
        
        Rscript scripts/run_cutadapt_CL.R \
            --input_file {input.reads_len_filter} \
            --primer_yaml {input.primer_yaml} \
            --output_path_linked "{output.linked}" \
            --output_path_unlinked "{output.unlinked}" \
            --cutadapt_error_rate "{CUTADAPT_ERROR}" \
            --n_cores 0

        seqkit rmdup "{output.linked}" > "{output.linked_dd}"

        chmod +x scripts/get_read_len_distribs.sh
        scripts/get_read_len_distribs.sh "{output.linked_dd}" "{output.linked_dd_distrib}" >> {log.trim_primers_log} 2>&1

        echo "Removing IDs of linked reads from unlinked reads" > {log.trim_primers_log}

        ## get IDs
        awk 'NR%4==1 {{sub(/^@/, "", $1); print $1}}' \
        "{output.linked_dd}" > "{output.linked_dd}.ids.txt"

        ## rm ids of linked from unlinked & deduplicate
        seqkit grep --invert-match -f "{output.linked_dd}.ids.txt" \
        "{output.unlinked}" | \
        seqkit rmdup > "{output.unlinked_sub}"
        """
