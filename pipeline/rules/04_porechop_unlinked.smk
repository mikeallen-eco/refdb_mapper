# rules/04_porechop_unlinked.smk

QUALITY_SCORE = config.get("quality_score", 10)
MIN_LENGTH = config.get("raw_min_length", 100)
MAX_LENGTH = config.get("raw_max_length", 500)
OUTPUT_DIR = config["paths"]["output"]

rule porechop_unlinked:
    input:
        reads_unlinked_sub = f"{OUTPUT_DIR}/03_trim_primers/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.unlinked.sub.fastq"
    output:
        reads_unlinked_sub_shuffle = temp(f"{OUTPUT_DIR}/04_porechop_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.unlinked.sub.shuf.fastq"),
        reads_porechopped = temp(f"{OUTPUT_DIR}/04_porechop_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.unlinked.sub.shuf.porechop.fastq"),
        reads_unlinked_distrib = temp(f"{OUTPUT_DIR}/04_porechop_unlinked/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.unlinked.sub.shuf.porechop.distrib.txt")
    conda:
        "envs/environment.yaml"
    log:
        porechop_unlinked_log = f"{OUTPUT_DIR}/04_porechop_unlinked/logs/{{sample}}_porechop_unlinked.log"
    shell:
        r"""
        echo "Porechopping unlinked reads"

        mkdir -p "{OUTPUT_DIR}/04_porechop_unlinked"

        # shuffle reads
        shuffle.sh in="{input.reads_unlinked_sub}" out="{output.reads_unlinked_sub_shuffle}"

        # run porechop with --discard_middle flag
        porechop_abi --threads 1 --discard_middle -i "{output.reads_unlinked_sub_shuffle}" > "{output.reads_porechopped}"

        # get read length distribution
        chmod +x scripts/get_read_len_distribs.sh
        scripts/get_read_len_distribs.sh "{output.reads_porechopped}" "{output.reads_unlinked_distrib}" >> {log.porechop_unlinked_log} 2>&1

        """
