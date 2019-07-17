# Thalassaemia point mutation/small InDel detection

This workflow describes: **(1)** how to rescue the poorly aligned NGS reads in Bam files, and **(2)** workflow for thalassaemia point mutation/small InDel detection on clusters managed by Portable Batch System(PBS).

First of all, if you start with large files such as whole genome sequencing BAM, you are suggested to only extract hemoglobin regions from BAM files by running the following commands:

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

And Python module dependency: pysam, numpy
____________________________________________________________________________________________________________

Now let's start, and suppose our working directory is /home/data/Thalaproject.

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
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step3_v1_HardFiltering_model.pbs
    |    |    |    |-- Poor_reads_rescue
    |    |    |    |    | -- Thalassemia.py
    |    |    |    |    | -- BamOPR.py
    |    |    |    |    | -- README.md
    |    |    |    |-- Known_Causal_Mutation
    |    |    |    |    | -- substitution.query_vars3.bed
    |    |    |    |    | -- indel.query_vars3.bed

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
    |    |    |    |    | -- Thala_Rescue_PBS
    |    |    |    |    |    | -- Thala_Rescue_Bam_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step2_GTing_model.pbs
    |    |    |    |    |    | -- Thala_Rescue_phase2_Step3_v1_HardFiltering_model.pbs
    |    |    | -- Bam_file
    |    |    |    | -- Sample1
    |    |    |    |    | -- Rescue_phase_Sample1.pbs
    |    |    |    | -- Sample2
    |    |    |    |    | -- Rescue_phase_Sample2.pbs
    |    |    |    | -- Sample3
    |    |    |    |    | -- Rescue_phase_Sample3.pbs
    |    |    | -- VCF_file
    |    |    |    | -- Sample1
    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_Sample1.pbs
    |    |    |    | -- Sample2
    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_Sample2.pbs
    |    |    |    | -- Sample3
    |    |    |    |    | -- Thala_Rescue_phase2_Step1_RunHC_Sample3.pbs
    |    |    |    | -- Joint
    |    |    | -- Thala_Rescue_phase2_Step2_GTing.pbs
    |    | -- submit.sh

* Step5: submit PBS files step-by-step by changing and running the following commands:

      vi submit.sh(change the target PBS scripts)
      sh submit.sh

This workflow is designed for clusters managed by PBS, for non-PBS servers, users can still run these scripts in bash(sh):

    cd Bam_file/Sample1
    sh Rescue_phase_Sample1.pbs
