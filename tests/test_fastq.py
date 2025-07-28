from src.kml_panelqc.fastq import prepare_fastq_by_samptab
from pathlib import Path

# input_tab = '/data/mengxf/Project/KML250721-cyp2d6-igtk/work/input-table/igtk-samples.tsv'
input_tab = '/data/mengxf/GitHub/KML-PanelQC/template/test.tsv'
work_dir = Path('/data/mengxf/Project/KML250721-cyp2d6-igtk/results/250725')

prepare_fastq_by_samptab(work_dir, input_tab, 32)
