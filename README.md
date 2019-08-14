# Thalassaemia point mutation/small InDel detection

This workflow describes: **(1)** how to rescue the poorly aligned NGS reads in Bam files, and **(2)** how to pick out thalassaemia point mutation/small InDel detection on clusters managed by Portable Batch System(PBS).

For PBS-free servers, users can still run these scripts in bash(sh) like:

    cd Bam_file/Sample1
    sh Rescue_phase_Sample1.pbs


First of all, if you start with large files such as whole genome sequencing BAM, you are strongly suggested to only extract hemoglobin regions from BAM files by running the following commands:

      samtools view -h -L Thalassaemia_hg19_genome.bed -b -o output.bam input.bam

## Prerequisite:

* Step0: Order your raw Bam files by **sample names** in the following strucure:

        | -- Bam_file_folder
        |   | -- Sample1
        |   |   | -- Sample1.bam
        |   | -- Thala_2
        |   |   | -- Thala_2.bam
        |   | -- TJLE
        |   |   | -- TJLE.bam

Folder names and Bam file names should be same, which should be the sample names.

Python module dependency: pysam, numpy, and commonly used software and their resources such as [GATK](https://software.broadinstitute.org/gatk/), [GATK_Bundle](https://software.broadinstitute.org/gatk/download/bundle), [Picard](https://broadinstitute.github.io/picard/) and [vcftools](https://vcftools.github.io/examples.html)
____________________________________________________________________________________________________________

Now let's start.

Suppose our working directory is **/home/data/Thalaproject**.

* Step1: create the rescue folder by running the following command in bash:

      cd /home/data/Thalaproject
      mkdir Rescue_Phase

* Step2: download this repository to the Rescue_Phase folder by running the following commands:

      cd /home/data/Thalaproject/Rescue_Phase
      git clone https://github.com/JavenCao/Thala_Rescue_workflow.git

now you should have the follwing structure:

    | -- /home/data/Thalaproject/
    |    | -- Rescue_Phase
    |    |    | -- Thala_Rescue_workflow
    |    |    |    | -- Thala_rescue_PBS.py
    |    |    |    | -- Thala_rescue_configuration.txt
    |    |    |    | -- supportingFun.py
    |    |    |    |    | -- PBSModels
    |    |    |    |    |    | -- Thala_Rescue_Bam_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step2_GTing_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step3_HardFiltering_model.pbs
    |    |    |    |    |    | -- Select_Annotation.pbs
    |    |    |    |    |    | -- plot.R    
    |    |    |    |    |    | -- Thala_Find_Causal_model.pbs
    |    |    |    |    |    | -- Find_causal.py
    |    |    |    |-- Poor_reads_rescue
    |    |    |    |    | -- Thalassemia.py
    |    |    |    |    | -- BamOPR.py
    |    |    |    |    | -- README.md
    |    |    |    |-- Known_Causal_Mutation
    |    |    |    |    | -- Known_InDel.vcf
    |    |    |    |    | -- Known_substitution.vcf

* Step3: go into the Thala_Rescue_workflow folder, and set parameters in the follwing file. The parameters are self-explainable.

      vi Thala_rescue_configuration.txt

* Step4: run the following commands in bash:

      python Thala_rescue_PBS.py

And after that you should have the following structure:

    | -- /home/data/Thalaproject/
    |    | -- Rescue_Phase
    |    |    | -- Thala_Rescue_workflow
    |    |    |    | -- Thala_rescue_PBS.py
    |    |    |    | -- Thala_rescue_configuration.txt
    |    |    |    | -- supportingFun.py
    |    |    |    |    | -- PBSModels
    |    |    |    |    |    | -- Thala_Rescue_Bam_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step2_GTing_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step3_HardFiltering_model.pbs
    |    |    |    |    |    | -- Select_Annotation.pbs
    |    |    |    |    |    | -- plot.R    
    |    |    |    |    |    | -- Thala_Find_Causal_model.pbs
    |    |    |    |    |    | -- Find_causal.py
    |    |    |    |-- Poor_reads_rescue
    |    |    |    |    | -- Thalassemia.py
    |    |    |    |    | -- BamOPR.py
    |    |    |    |    | -- README.md
    |    |    |    |-- Known_Causal_Mutation
    |    |    |    |    | -- Known_InDel.vcf
    |    |    |    |    | -- Known_substitution.vcf
    --------------------------------------------------------------------------
    |    |    | -- Bam_file                                                  |
    |    |    |    | -- Sample1                                              |
    |    |    |    |    | -- Rescue_phase_Sample1.pbs                        |
    |    |    |    | -- Sample2                                              |
    |    |    |    |    | -- Rescue_phase_Sample2.pbs                        |
    |    |    |    | -- Sample3                                              |
    |    |    |    |    | -- Rescue_phase_Sample3.pbs                        |
    |    |    | -- VCF_file                                                  |
    |    |    |    | -- Sample1                                              |
    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_Sample1.pbs     |
    |    |    |    | -- Sample2                                              |
    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_Sample2.pbs     |
    |    |    |    | -- Sample3                                              |
    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_Sample3.pbs     |
    |    |    |    | -- Joint                                                |
    |    |    |    |    | -- Thala_Rescue_phase2_Step3_HardFiltering.pbs     |
    |    |    |    |    | -- Select_Annotation.pbs                           |
    |    |    |    |    | -- plot.R                                          |
    |    |    |    |    | -- Find_Causal.pbs                                 |
    |    |    |    |    | -- Find_causal.py                                  |
    |    |    | -- Thala_Rescue_phase2_Step2_GTing.pbs                       |
    |    | -- submit.sh                                                      |
    --------------------------------------------------------------------------
* Step5: submit PBS files step-by-step by changing and running the following commands:

      vi submit.sh(step-by-step, change the target PBS scripts)
      sh submit.sh

After get gVCF for each sample, go to the VCF_file folder, and do Joint Genotyping and hard filtering by running the following command:

    cd ./VCF_file
    qsub Thala_Rescue_phase2_Step2_GTing.pbs
    cd ./Joint
    qsub Thala_Rescue_phase2_Step3_HardFiltering.pbs

For the **Hard filtering steps**, we suggest to set cut-offs based on the emprical distribution of each annotation.
Here, we provide scripts to select and plot the emprical distribution for each annotation.

    Select_Annotation.pbs
    plot.R

To increase sensitivity, the default cut-offs are as follows:

      SNP:"QD<2.0 || FS>60 || MQ<40 || MQRankSum<-12.5 || ReadPosRankSum<-8.0"
      INDEL:"QD < 2.0 || FS > 200 || ReadPosRankSum < -20 || SOR > 10.0"

After Step5, in Joint folder, you should get SNP and INDEL vcf files for further process.

* Step6: Find current known thalassaemia causal mutations based on [HbVar](http://globin.cse.psu.edu/hbvar/menu.html) and [ITHANET](https://www.ithanet.eu/). We have created collections from these databases, and you just running the follwing commands to pick these mutations out.

      cd ./VCF_file/Joint/
      qsub Find_Causal.pbs (or sh Find_causal.pbs for non-PBS servers)

After this, you will find a newly generated folder named **Prediction**, only Positive samples are reported with mutation details.

This workflow is designed for clusters managed by PBS, for PBS-free servers, users can still run these scripts in bash(sh) like:

    cd Bam_file/Sample1
    sh Rescue_phase_Sample1.pbs

## License

This project is licensed under GNU GPL v3.

## Authors

Cao Yujie(The University of Hong Kong)


