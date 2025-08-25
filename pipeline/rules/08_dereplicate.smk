# rules/08_dereplicate.smk

QUALITY_SCORE = config.get("quality_score", 10)
MIN_LENGTH = config.get("raw_min_length", 100)
MAX_LENGTH = config.get("raw_max_length", 500)
MIN_LENGTH_FIN = config.get("fin_min_length", 100)
MAX_LENGTH_FIN = config.get("fin_max_length", 500)
OUTPUT_DIR = config["paths"]["output"]

rule dereplicate:
    """
    Dereplicate reads per sample (exact matches) and output FASTA with size counts.
    Sequence IDs are prefixed with sample name.
    """
    input:
        reads_fasta_att = f"{OUTPUT_DIR}/07_simplify_fasta/{{sample}}.q{QUALITY_SCORE}.l{MIN_LENGTH}.L{MAX_LENGTH}.linked_unlinked.sub.dd.l{MIN_LENGTH_FIN}.L{MAX_LENGTH_FIN}.att.fasta",
    output:
        derep_fasta = temp(f"{OUTPUT_DIR}/08_dereplicate/{{sample}}_dereplicate.fasta")
    conda:
        "envs/environment.yaml"
    log:
        f"{OUTPUT_DIR}/08_dereplicate/logs/{{sample}}_dereplicate.log"
    shell:
        r"""
        mkdir -p {OUTPUT_DIR}/08_dereplicate

        vsearch \
            --derep_fulllength {input.reads_fasta_att} \
            --output {output.derep_fasta} \
            --sizeout \
            --relabel {wildcards.sample}_ \
            > {log} 2>&1

        readcount=$(grep -o 'size=[0-9]*' {output.derep_fasta} \
            | cut -d= -f2 \
            | paste -sd+ - \
            | bc)
        echo "File contains ${{readcount}} total reads (accounting for ;size=)" >> {log}

        """
