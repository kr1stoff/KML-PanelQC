rule copy_reference:
    input:
        ref=config["reference"],
    output:
        temp("prepare/ref.fa"),
    benchmark:
        ".log/prepare/copy_reference.bm"
    log:
        ".log/prepare/copy_reference.log",
    shell:
        "cp {input.ref} {output}"


rule prepare_bwa_index:
    input:
        rules.copy_reference.output,
    output:
        multiext(rules.copy_reference.output[0], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    benchmark:
        ".log/prepare/prepare_bwa_index.bm"
    log:
        ".log/prepare/prepare_bwa_index.log",
    conda:
        config["conda"]["basic"]
    shell:
        "bwa index {input}"


rule bedtools_sort:
    output:
        "prepare/target.sorted.bed",
    benchmark:
        ".log/prepare/bedtools_sort.bm"
    log:
        ".log/prepare/bedtools_sort.log",
    conda:
        config["conda"]["basic2"]
    shell:
        "bedtools sort -i {config[bed]} > {output} 2> {log}"
