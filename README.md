# KML-PanelQC

试剂盒性能评估, 主要指标为 数据量/Q30/重复率/比对率/在靶率/50X覆盖度/均一性 等

## 设计思路

1. fastp 结果统计 fastq 质量
2. samtools stat 统计比对率, `-t` 参数指定靶区域获取在靶率
3. samtools depth 指定 `bed` 区域获取 4X/10X/30X/50X/100X 的覆盖率, 以及 20%/50% 均一性
4. 比对率/在靶率使用 read 数量计算, 覆盖度/覆盖深度/均一性使用 base 数量计算

## 命令行

- 主程序
  1.浩博 HBV

  ```bash
  ~/miniforge3/envs/python3.12/bin/python -m src.kml_panelqc \
    --input-tab /data/mengxf/Project/KML250721-cyp2d6-igtk/work/input-table/igtk-samples.tsv \
    --output-dir /data/mengxf/Project/KML250721-cyp2d6-igtk/results/250725 \
    --reference /data/mengxf/Database/reference/hg38/hg38.fa \
    --bed /data/mengxf/Project/KML250721-cyp2d6-igtk/probeCov.gene.bed \
    --threads 32
  ```

- snakemake 运行

  ```bash
  mamba run -n snakemake \
    snakemake -c 32 --use-conda \
    -s /data/mengxf/GitHub/KML-PanelQC/wf-panelqc/Snakefile \
    --configfile .temp/snakemake.yaml \
    --scheduler greedy \
    --ignore-incomplete
  ```
