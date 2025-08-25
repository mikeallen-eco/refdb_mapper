# rules/05_combine_linked_unlinked.smk

QUALITY_SCORE = config.get("quality_score", 10)
MIN_LENGTH = config.get("raw_min_length", 100)
MAX_LENGTH = config.get("raw_max_length", 500)
OUTPUT_DIR = config["paths"]["output"]

rule combine_linked_unlinked:
    input:
        reads_linked_dd = f"{OUTPUT_DIR}/03_trim_primers/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked.dd.fastq",
        reads_porechopped = f"{OUTPUT_DIR}/04_porechop_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.unlinked.sub.shuf.porechop.fastq"
    output:
        reads_linked_unlinked = temp(f"{OUTPUT_DIR}/05_combine_linked_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.fastq"),
        reads_linked_unlinked_dd = temp(f"{OUTPUT_DIR}/05_combine_linked_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.fastq"),
        reads_linked_unlinked_dd_distrib = f"{OUTPUT_DIR}/05_combine_linked_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.distrib.txt"
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/05_combine_linked_unlinked/logs/{{sample}}_combine_linked_unlinked.log"
    shell:
        r"""
        echo "Combining linked and unlinked reads"

        mkdir -p "{OUTPUT_DIR}/05_combine_linked_unlinked"

        # Combine linked and porechopped unlinked reads
        cat "{input.reads_linked_dd}" "{input.reads_porechopped}" > "{output.reads_linked_unlinked}"

        # Remove duplicates
        seqkit rmdup "{output.reads_linked_unlinked}" -o "{output.reads_linked_unlinked_dd}"

        # Get read length distribution
        chmod +x scripts/get_read_len_distribs.sh
        scripts/get_read_len_distribs.sh "{output.reads_linked_unlinked_dd}" "{output.reads_linked_unlinked_dd_distrib}" >> {log} 2>&1
        """
