from subprocess import run
import pandas as pd
from typing import Generator
import logging


def check_ref_and_bed(reference: str, bed: str) -> None:
    """
    检查参考基因组和 BED 文件
    :reference: 参考基因组路径
    :bed: BED 文件路径
    """
    logging.info("检查参考基因组和 BED 文件")
    # reference id 列表
    ref_ids = list(get_ref_ids(reference))
    # bed 第一列
    bed_ids = pd.read_table(bed, header=None, usecols=[
                            0]).iloc[:, 0].unique().tolist()
    not_in_ref_ids = [bed_id for bed_id in bed_ids if bed_id not in ref_ids]
    assert len(
        not_in_ref_ids) == 0, f"BED 文件中存在不在参考基因组中的 ID: {', '.join(not_in_ref_ids)}"


def get_ref_ids(reference: str) -> Generator:
    """
    获取参考基因组 ID 列表
    :reference: 参考基因组路径
    :return: 参考基因组 ID 生成器
    """
    with open(reference) as f:
        for line in f:
            if line.startswith('>'):
                ref_id = line.strip('>').split()[0]
                yield ref_id
