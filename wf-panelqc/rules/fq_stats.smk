rule fastp_pe:
    input:
        ".rawdata/{sample}_1.fastq.gz",
        ".rawdata/{sample}_2.fastq.gz",
    output:
        j="fq_stats/{sample}.json",
        h="fq_stats/{sample}.html",
        o=temp("fq_stats/{sample}.1.fastq.gz"),
        O=temp("fq_stats/{sample}.2.fastq.gz"),
    benchmark:
        ".log/fq_stats/{sample}.fastp_pe.bm"
    log:
        ".log/fq_stats/{sample}.fastp_pe.log",
    conda:
        config["conda"]["basic2"]
    threads: config["threads"]["low"]
    shell:
        "fastp -w {threads} -j {output.j} -h {output.h} -o {output.o} -O {output.O} -i {input[0]} -I {input[1]} &> {log}"


rule fq_stats_summary:
    input:
        expand("fq_stats/{sample}.json", sample=config["samples"]),
    output:
        "fq_stats/fq_summary.tsv",
    benchmark:
        ".log/fq_stats/fq_all_samples_qc.bm"
    log:
        ".log/fq_stats/fq_all_samples_qc.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/fq_all_samples_qc.py"
