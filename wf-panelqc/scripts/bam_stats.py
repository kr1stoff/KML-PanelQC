import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], 'w')


def get_bam_stat_dict(bam_stat_file) -> dict:
    bam_stat_dict = {}
    with open(bam_stat_file, 'r') as f:
        for line in f:
            lns = line.split('#')[0].strip().split(':\t')
            bam_stat_dict[lns[0]] = lns[1]
    return bam_stat_dict


# 靶区域 bam stat 字典
target_stats = get_bam_stat_dict(snakemake.input[1])
# 全部比对结果 bam stat 字典
all_stats = get_bam_stat_dict(snakemake.input[0])
# 比对率 reads mapped and paired / raw total sequences
mapped_reads = int(all_stats['reads mapped and paired'])
mapped_rate = mapped_reads / int(all_stats['raw total sequences'])
# 在靶率 (target) reads mapped and paired / (all) reads mapped and paired
ontarget_reads = int(target_stats['reads mapped and paired'])
ontarget_rate = ontarget_reads / mapped_reads

target_depth = pd.read_table(snakemake.input[2], header=None)
# 覆盖率
target_size = target_depth.shape[0]
target_covered_size = target_depth[target_depth[2] > 0].shape[0]
covarage_rate = target_covered_size / target_size
# target_size, target_covered_size, covarage_rate

# 4X/10X/30X/50X/100X 覆盖深度
cover4x_rate = target_depth[target_depth[2] > 4].shape[0] / target_size
cover10x_rate = target_depth[target_depth[2] > 10].shape[0] / target_size
cover30x_rate = target_depth[target_depth[2] > 30].shape[0] / target_size
cover50x_rate = target_depth[target_depth[2] > 50].shape[0] / target_size
cover100x_rate = target_depth[target_depth[2] > 100].shape[0] / target_size
# 20%/50% 均一性
mean_depth = int(target_depth[2].mean())
depth_20_rate = target_depth[target_depth[2] > mean_depth * 0.2].shape[0] / target_size
depth_50_rate = target_depth[target_depth[2] > mean_depth * 0.5].shape[0] / target_size

# * 输出
outputs = [mapped_reads, mapped_rate, ontarget_reads, ontarget_rate, target_size, target_covered_size, covarage_rate,
           cover4x_rate, cover10x_rate, cover30x_rate, cover50x_rate, cover100x_rate, mean_depth, depth_20_rate, depth_50_rate]
outputs = [o if type(o) == int else round(o, 4) for o in outputs]
df = pd.DataFrame([outputs], columns=['MappedReads', 'MappedRate', 'OnTargetReads', 'OnTargetRate', 'TargetSize', 'TargetCoveredSize', 'CoverageRate',
                  '4xCoverageRate', '10xCoverageRate', '30xCoverageRate', '50xCoverageRate', '100xCoverageRate', 'MeanDepth', '20xDepthRate', '50xDepthRate'])
df.to_csv(snakemake.output[0], index=False)
