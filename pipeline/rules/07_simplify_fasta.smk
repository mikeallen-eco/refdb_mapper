# rules/07_simplify_fasta.smk

QUALITY_SCORE = config.get("quality_score", 10)
MIN_LENGTH = config.get("raw_min_length", 100)
MAX_LENGTH = config.get("raw_max_length", 500)
MIN_LENGTH_FIN = config.get("fin_min_length", 100)
MAX_LENGTH_FIN = config.get("fin_max_length", 500)
OUTPUT_DIR = config["paths"]["output"]

rule simplify_fasta:
    input:
        reads_final_length_filter = f"{OUTPUT_DIR}/06_final_length_filter/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.l{MIN_LENGTH_FIN}.L{MAX_LENGTH_FIN}.fastq",
    output:
        reads_final_length_filter_fasta = temp(f"{OUTPUT_DIR}/07_simplify_fasta/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.l{MIN_LENGTH_FIN}.L{MAX_LENGTH_FIN}.fasta"),
        reads_final_length_filter_fasta_att = temp(f"{OUTPUT_DIR}/07_simplify_fasta/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.l{MIN_LENGTH_FIN}.L{MAX_LENGTH_FIN}.att.fasta"),
    conda:
        "envs/environment.yaml"
    log:
        simplify_fasta_log = f"{OUTPUT_DIR}/07_simplify_fasta/logs/{{sample}}_simplify_fasta.log"

    shell:
        r"""
        echo "Creating fasta file with a minimal header." >> {log.simplify_fasta_log} 2>&1

        mkdir -p "{OUTPUT_DIR}/07_simplify_fasta"

        seqkit fq2fa "{input.reads_final_length_filter}" -o "{output.reads_final_length_filter_fasta}"

        ## Drop attributes except barcode
        awk 'BEGIN {{OFS=" "}} /^>/ {{printf "%s", $1; for(i=2; i<=NF; i++) if($i ~ /^barcode=/) printf " %s", $i; print ""}} !/^>/ {{print}}' "{output.reads_final_length_filter_fasta}" > "{output.reads_final_length_filter_fasta_att}"

        # count reads and save to log file
        seqkit stats "{output.reads_final_length_filter_fasta_att}" >> {log.simplify_fasta_log} 2>&1
        """
