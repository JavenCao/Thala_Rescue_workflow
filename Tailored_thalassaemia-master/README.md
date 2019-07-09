# Tailored_thalassaemia

This module(including Thalassemia.py and BamOPR.py) was written to rescue the poorly aligned NGS reads falling into the homologous hemoglobin gene clusters, such as HBA1 and HBA2. In short, we first define a list of poorly aligned reads, then by random assignment of these paired-end reads while keeping reasonable insert size distribution, we rescued the raw Bam file and generated modified Bam file for thalassaemia point mutation/small InDel detection.

This module can be used to rescue any Bam files, such as WGS-bam, WES-bam, targeted-seq-Bam.

A quick look of the usage is demonstrated below:

    usage: Thalassemia.py [-h] [-bam BAMFILE] [-o OUTPUT] [-ref REFERENCE]

    optional arguments:
      -h, --help            show this help message and exit
      -bam BAMFILE, --bamfile BAMFILE
                            /Path/to/input_Bam_file.bam, example:/home/data/input.bam
      -o OUTPUT, --output OUTPUT
                            /Path/to/output_Bam_file.bam,example:/home/data/out.bam
      -ref REFERENCE, --reference REFERENCE
                            reference of human genome assebly, hg19(Default) or GRCh38

Notes: Python module dependency:

    pysam, random, numpy

Also, a typical workflow from rescue Bam to variants calling based on GATK is illustrated [here](https://github.com/JavenCao/Thala_Rescue_workflow).
