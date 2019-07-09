# Thala_Rescue_workflow

This workflow describes how to rescue the Bam files for thalassaemia mutation detecton in a batch mode on clusters managed by PBS.

Pre-request:

Run workflow of Easy_WES_By_PBS(at least get bam files in each sample folder).

Now let's start the rescue process, and suppose our working directory is /home/data/Thalaproject.

Step1: create the rescue folder by running the following command in bash:

    cd /home/data/Thalaproject
    mkdir Rescue_Phase

Step2: upload this folder to /home/data/WESprojec/Rescue_Phase, or run the following commands:

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
    |    |    | -- Rescue_code
    |    |    |    | -- file1.pbs
    |    |    |    | -- file2.pbs
    |    |    |    | -- file3.pbs

Step3: go into the Rescue_code, and set parameters for the rescue process.

Step4: run the following commands in bash:

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
    |    |    | -- Rescue_code
    |    |    |    | -- file1.pbs
    |    |    |    | -- file2.pbs
    |    |    |    | -- file3.pbs
    |    |    | -- Bam_file
    |    |    |    | -- Sample1
    |    |    |    |    | -- Sample1.rescue.pbs
    |    |    |    | -- Sample2
    |    |    |    |    | -- Sample2.rescue.pbs
    |    |    |    | -- Sample3
    |    |    |    |    | -- Sample3.rescue.pbs
    |    |    | -- VCF_file
    |    |    |    | -- Sample1
    |    |    |    |    | -- Thala_Rescue_RunHC_Sample1.pbs
    |    |    |    | -- Sample2
    |    |    |    |    | -- Thala_Rescue_RunHC_Sample2.pbs
    |    |    |    | -- Sample3
    |    |    |    |    | -- Thala_Rescue_RunHC_Sample3.pbs
    |    |    |    | -- Joint
    |    |    | -- Easy_WES_jointGT.pbs
    |    | -- submit.sh

Step5: submit PBS files step-by-step by changing and running the following commands:

    vi submit.sh(change the target PBS scripts)
    sh submit.sh
