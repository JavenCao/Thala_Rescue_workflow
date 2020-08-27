import os
from supportingFun import *

Path_dict = {}
SampleList = []

# load values from config file
Path_dict, SampleList = load_config_file("Thala_rescue_configuration.txt")

# Create folders
wkd = Path_dict.get('rescue_folder')
#Raw_data_folder = wkd+'/Raw_data'
Bam_file_folder = wkd + '/Bam_file'
VCF_file_folder = wkd + '/VCF_file'
Create_Folders(Bam_file_folder, SampleList)
Create_Folders(VCF_file_folder, SampleList)

#-------------------------------------begin GATK process----------------------------------------
SampleNum = len(SampleList)
PBSModel_folder = os.getcwd() + '/PBSModels'

# Create rescue PBS files in each sample folder
rescue_model_file = PBSModel_folder + '/' + Path_dict.get('PBSfile1')
ModifyAndCreate(rescue_model_file, Path_dict,
                Bam_file_folder, SampleList, "Rescue_phase")

# For running HC and getting RawVariants.g.vcf files
phase2_RunHC_model_file = PBSModel_folder + '/' + Path_dict.get('PBSfile2')
ModifyAndCreate(phase2_RunHC_model_file, Path_dict,
                VCF_file_folder, SampleList, "Thala_Rescue_phase2_Step1_RunHC")

# For Jointing Genotyping
Joint_folder = VCF_file_folder + '/Joint'
os.mkdir(Joint_folder)
# copy plot.R for annotation plotting purpose
plot_R_file = PBSModel_folder + '/plot.R'
command_line = """cp {0} {1}""".format(plot_R_file, Joint_folder)
os.system(command_line)
# copy Find_causal.py for finding causal mutation purpose
Find_causal_file = PBSModel_folder + '/Find_causal.py'
command_line = """cp {0} {1}""".format(Find_causal_file, Joint_folder)
os.system(command_line)

phase2_GTing_model_file = PBSModel_folder + '/' + Path_dict.get('PBSfile3')
ModifyAndCreate_v2(phase2_GTing_model_file, Path_dict,
                   VCF_file_folder, SampleList, 'Thala_Rescue_phase2_Step2_GTing')

# Hard filtering
phase2_HardF_model_file = PBSModel_folder + '/' + Path_dict.get('PBSfile4')
ModifyAndCreate_v2(phase2_HardF_model_file, Path_dict, Joint_folder,
                   SampleList, 'Thala_Rescue_phase2_Step3_HardFiltering')

#-------------------------------------end GATK process----------------------------------------


# create submit.sh file under project folder for submit tasks to CGS clusers
submitFile = wkd + '/submit.sh'
t = ''
for sample in SampleList:
  t = t + sample + ' '
line1 = "for i in " + t + '\n'
otherlines = "do\n" + "\tcd " + wkd + "/Bam_file/$i" + '\n' + \
    "\tqsub Rescue_phase_" + "$i" + ".pbs\n" + \
    "#\tcd " + wkd + "/VCF_file/$i" + '\n' + \
    "#\tqsub Thala_Rescue_phase2_Step1_RunHC_" + "$i" + ".pbs\n" + "done\n"

with open(submitFile, "w") as subF:
  subF.write(line1)
  subF.write(otherlines)
