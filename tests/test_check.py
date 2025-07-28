from src.kml_panelqc.check import check_ref_and_bed


bed = '/data/mengxf/GitHub/KML-CYP2D6/bed/probeCov.gene.bed'
# bed = '/data/mengxf/GitHub/KML-PanelQC/tests/probeCov.gene.bed'
reference = '/data/mengxf/Database/reference/hg38/hg38.fa'
check_ref_and_bed(reference, bed)
