# mbd56_pollen

Contains code used for Ichino et al. 2022 (Cell Reports, link: TBD) for updating gene annotations based on pollen RNA-seq. We primarily used the Mikado pipeline (https://mikado.readthedocs.io/en/stable/), but include here an additional script that was used to refine the predictions selected by Mikado to better detect the non-coding sequences and reactivated TEs we were seeing in this tissue.

Contains 3 files:
- mikado_refine.py - this script refined the output of Mikado; instructions for how to run included in file header, including a detailed overview of how to include this in a 'regular' Mikado run
- mikado_refine_version2_final_wCM.gtf - the final updated annotations, including all transcripts passing filters
- mikado_refine_version2_final_noAS_wCM.gtf - final updated annotations, but with annotations that were noncoding + antisense to protein-coding genes removed

Script dependencies (python v.3.x):
- argparse
- pysam
- pyBigWig
- numpy
- Bio 
