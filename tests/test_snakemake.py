from src.kml_panelqc.snakemake import create_snakemake_configfile
from pathlib import Path

# input_tab = '/data/mengxf/Project/KML250721-cyp2d6-igtk/work/input-table/igtk-samples.tsv'
input_tab = '/data/mengxf/GitHub/KML-PanelQC/template/test.tsv'
work_dir = Path('/data/mengxf/Project/KML250721-cyp2d6-igtk/results/250725')
bed = '/data/mengxf/GitHub/KML-CYP2D6/bed/probeCov.gene.bed'
reference = '/data/mengxf/Database/reference/hg38/hg38.fa'

create_snakemake_configfile(input_tab, work_dir, 32, reference, bed)
