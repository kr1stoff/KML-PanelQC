rule bwa_mem:
    input:
        ".rawdata/{sample}_1.fastq.gz",
        ".rawdata/{sample}_2.fastq.gz",
    output:
        "bam_stats/{sample}.bam",
    benchmark:
        ".log/bam_stats/{sample}.bwa_mem.bm"
    log:
        ".log/bam_stats/{sample}.bwa_mem.log",
    conda:
        config["conda"]["basic"]
    params:
        bwa="-M -Y -R '@RG\\tID:{sample}\\tSM:{sample}'",
        view="-hbS",
    threads: config["threads"]["high"]
    shell:
        """
        bwa mem -t {threads} {params.bwa} {config[reference]} {input} 2> {log} | \
            samtools view -@ {threads} {params.view} - 2>> {log} | \
            samtools sort -@ {threads} -o {output} - 2>> {log}
        """


rule samtools_stat:
    input:
        rules.bwa_mem.output,
        ".temp/target.sorted.bed",
    output:
        all_stat="bam_stats/{sample}.bam.stat",
        target_stat="bam_stats/{sample}.bam.target.stat",
    benchmark:
        ".log/bam_stats/{sample}.samtools_stat.bm"
    log:
        ".log/bam_stats/{sample}.samtools_stat.log",
    conda:
        config["conda"]["basic"]
    threads: config["threads"]["low"]
    shell:
        """
        samtools stat {input[0]} | grep ^SN | cut -f 2- > {output.all_stat} 2> {log}
        samtools stat -t {input[1]} {input[0]} | grep ^SN | cut -f 2- > {output.target_stat} 2>> {log}
        """


rule samtools_depth:
    input:
        rules.bwa_mem.output,
        ".temp/target.sorted.bed",
    output:
        "bam_stats/{sample}.bam.target.depth",
    benchmark:
        ".log/bam_stats/{sample}.samtools_depth.bm"
    log:
        ".log/bam_stats/{sample}.samtools_depth.log",
    conda:
        config["conda"]["basic"]
    shell:
        "samtools depth -b {input[1]} -a {input[0]} -o {output} 2> {log}"


rule bam_stats:
    input:
        rules.samtools_stat.output.all_stat,
        rules.samtools_stat.output.target_stat,
        rules.samtools_depth.output,
    output:
        "bam_stats/{sample}.stats.csv",
    benchmark:
        ".log/bam_stats/{sample}.bam_stats.bm"
    log:
        ".log/bam_stats/{sample}.bam_stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats.py"


rule bam_stats_summary:
    input:
        expand("bam_stats/{sample}.stats.csv", sample=config["samples"]),
    output:
        "bam_stats/bam_summary.tsv",
    benchmark:
        ".log/bam_stats/bam_stats_summary.bm"
    log:
        ".log/bam_stats/bam_stats_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats_summary.py"
