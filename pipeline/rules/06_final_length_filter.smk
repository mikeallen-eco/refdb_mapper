# rules/06_final_length_filter.smk

QUALITY_SCORE = config.get("quality_score", 10)
MIN_LENGTH = config.get("raw_min_length", 100)
MAX_LENGTH = config.get("raw_max_length", 500)
MIN_LENGTH_FIN = config.get("fin_min_length", 100)
MAX_LENGTH_FIN = config.get("fin_max_length", 500)
OUTPUT_DIR = config["paths"]["output"]

rule final_length_filter:
    input:
        reads_linked_unlinked_dd = f"{OUTPUT_DIR}/05_combine_linked_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.fastq",
    output:
        final_reads_lenfilt = temp(f"{OUTPUT_DIR}/06_final_length_filter/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.l{MIN_LENGTH_FIN}.L{MAX_LENGTH_FIN}.fastq"),
        final_reads_lenfilt_distrib = f"{OUTPUT_DIR}/06_final_length_filter/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.l{MIN_LENGTH_FIN}.L{MAX_LENGTH_FIN}.distrib.txt"
    conda:
        "envs/environment.yaml"
    log:
        final_lengfilt_log = f"{OUTPUT_DIR}/06_final_length_filter/logs/{{sample}}_final_length_filter.log"
    shell:
        r"""
        echo "Performing final length filter: min {MIN_LENGTH_FIN}, max {MAX_LENGTH_FIN}."

        mkdir -p "{OUTPUT_DIR}/06_final_length_filter"

        seqkit seq --min-len {MIN_LENGTH_FIN} --max-len {MAX_LENGTH_FIN} {input.reads_linked_unlinked_dd} -o {output.final_reads_lenfilt} >> {log.final_lengfilt_log} 2>&1

        # Get read length distribution
        chmod +x scripts/get_read_len_distribs.sh
        scripts/get_read_len_distribs.sh "{output.final_reads_lenfilt}" "{output.final_reads_lenfilt_distrib}" >> {log.final_lengfilt_log} 2>&1
        """
