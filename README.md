# Thala_Rescue_workflow

This workflow describes how to rescue the Bam files for thalassaemia mutation detecion on clusters managed by PBS.

## Pre-request:

* Bam files in each sample folders, and the folder structure is illustrated [here](https://github.com/JavenCao/Easy_WES_By_PBS)(Step7)

* Python module dependency: pysam, numpy


Now let's start the rescue process, and suppose our working directory is /home/data/Thalaproject.

* Step1: create the rescue folder by running the following command in bash:

      cd /home/data/Thalaproject
      mkdir Rescue_Phase

* Step2: upload this folder to /home/data/WESprojec/Rescue_Phase, or run the following commands:

      cd /home/data/Thalaproject/Rescue_Phase
      git clone https://github.com/JavenCao/Thala_Rescue_workflow.git

now you should have the follwing structure:

    |    |/home/data/Thalaproject/
    |    | -- Raw_data
    |    | .... ignore here
    |    | -- Bam_file
    |    | .... ignore here
    |    | -- VCF_file
    |    | .... ignore here
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
    |    |    |    |-- Tailored_thalassaemia-master
    |    |    |    |    | -- Thalassemia.py
    |    |    |    |    | -- BamOPR.py
    |    |    |    |    | -- README.md
    |    |    |    |-- Known_Causal_Mutation
    |    |    |    |    | -- substitution.query_vars3.bed
    |    |    |    |    | -- indel.query_vars3.bed

* Step3: go into the Thala_Rescue_workflow folder, and set parameters in the follwing file. The parameters are self-explainable.

      Thala_rescue_configuration.txt

* Step4: run the following commands in bash:

      python rescue_thala.py

And after you should have the following structure:

    |    |/home/data/Thalaproject/
    |    | -- Raw_data
    |    | .... ignore here
    |    | -- Bam_file
    |    | .... ignore here
    |    | -- VCF_file
    |    | .... ignore here
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

This workflow is designed on clusters managed by PBS, for non-PBS servers, users are suggested to combine the PBS scritps into a single bash script before running the commands.
