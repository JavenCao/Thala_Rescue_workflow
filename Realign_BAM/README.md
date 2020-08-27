# Realignment of RMAs in homologous haemoglobin regions

This module (including Thalassemia.py and supproting functions - BamOPR.py) was disclosed to recover and realign reads with mulitple alignments(RMAs) in the homologous hemoglobin gene clusters, such as HBA2, HBA1 and their pseudogene, HBG1 and HBG2, for better detection of thalassaemia pathogenic point mutations and small InDels.

In short, we first define a list of RMAs that needs to be re-aligned, then by random assignment of these paired-end reads while keeping reasonable insert size, we rescued the raw Bam file and generated modified Bam file, which can be input for varaints callers, such as GATK-HaplotypeCaller.

This module can be used to rescue bwa aligned Bam files from whole genome sequencing(WGS), whole exome sequencing(WES) and panel-based targeted sequencing.

A quick look of the usage is shown below:

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

Example:

    python Thalassemia.py -h
    python Thalassemia.py -bam /home/data/input.bam -o /home/result/output.bam -ref "hg19"
## (A) Point mutation/small InDel

For point mutations and small InDels, a typical workflow starting from RMA realignment to variants calling, finally to identify thalassaemia pathogenic variants are illustrated [here](https://github.com/JavenCao/Thala_Rescue_workflow).

## (B) Structural variantion

For SV, workflow is illustrated [here](https://github.com/JavenCao/Tailored_SV_thala).

## License

This project is licensed under GNU GPL v3.

## Authors

Cao Yujie(The University of Hong Kong)
