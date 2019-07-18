import os


def load_config_file(config_name):
    """ load config values from config file"""
    #
    config_var = ['rescue_folder', 'Raw_Bam_file_folder', 'BWA_path', 'samtools_path', 'picard_path', 'GATK_path', 'GATK_bundle_path',
                  'ANNO_path', 'queue', 'walltime', 'nodes', 'ppn', 'mem', 'Email', 'PBSfile1', 'PBSfile2', 'PBSfile3', 'PBSfile4', 'PBSfile5']

    config_dict = {}

    SampleList = []

    with open(config_name) as ConfigFile:
        for line in ConfigFile:
            line = line.strip('\n')
            # Sample lines
            if(line.startswith("sample")):
                SampleInfo = line.split('=')
                SampleList.append(SampleInfo[1])
            else:
                info = line.split('=')
                if(info[0]) in config_var:
                    config_dict[info[0]] = info[1]

    return config_dict, SampleList


def Create_Folders(OuterSide, InnerListNames):
    """First create OuterSide folder, then Within the OuterSide, create each Inner folder in the order of InnerListNanes"""
    try:
        os.mkdir(OuterSide)
        for i in InnerListNames:
            subfolder = OuterSide + '/' + i
            os.mkdir(subfolder)
    except OSError:
        print """
        |---------------------------------------|
        | Please delete all the newly created   |
        | folders and re-run the program again !|
        |---------------------------------------|
        """
    return 1


def ModifyAndCreate(modelfile, Path_dict, Outer_folder, SampleList, prefix):
    """under each sample folder, create its own PBS files, with name prefix_sample.pbs"""
    with open(modelfile) as file:
        l = file.readlines()
        # strip the tailing '\n' of each l element
        # l = [ l[i].strip('\n') for i in range(len(l))]
        for sample in SampleList:

            for i in range(len(l)):

                if(l[i].startswith("#PBS -N")):
                    l[i] = "#PBS -N " + sample + '_' + prefix + '\n'

                elif(l[i].startswith("#PBS -l")):
                    l[i] = """#PBS -l mem={0},nodes={1}:ppn={2},walltime={3}{4}""".format(Path_dict.get(
                        'mem'), Path_dict.get('nodes'), Path_dict.get('ppn'), Path_dict.get('walltime'), '\n')

                elif(l[i].startswith("#PBS -q")):
                    l[i] = """#PBS -q {0}{1}""".format(
                        Path_dict.get('queue'), '\n')

                elif(l[i].startswith("#PBS -m abe -M")):
                    l[i] = "#PBS -m abe -M " + Path_dict.get('Email') + '\n'

                elif(l[i].startswith("i=")):
                    l[i] = "i=" + sample + '\n'

                elif(l[i].startswith("wkd=")):
                    l[i] = "wkd=" + Path_dict.get('rescue_folder') + '\n'

                elif(l[i].startswith("Raw_Bam_file_folder=")):
                    l[i] = "Raw_Bam_file_folder=" + \
                        Path_dict.get('Raw_Bam_file_folder') + '\n'

                elif(l[i].startswith("GATK_Bundle")):
                    l[i] = "GATK_Bundle=" + \
                        Path_dict.get("GATK_bundle_path") + '\n'

                elif(l[i].startswith("BWA")):
                    l[i] = "BWA=" + Path_dict.get("BWA_path") + '\n'

                elif(l[i].startswith("SAMTOOLS")):
                    l[i] = "SAMTOOLS=" + Path_dict.get("samtools_path") + '\n'

                elif(l[i].startswith("PICARD")):
                    l[i] = "PICARD=" + Path_dict.get("picard_path") + '\n'

                elif(l[i].startswith("GATK")):
                    l[i] = "GATK=" + Path_dict.get("GATK_path") + '\n'

                elif(l[i].startswith("ANNO")):
                    l[i] = "ANNO=" + Path_dict.get('ANNO_path')

            newFileName = Outer_folder + '/' + sample + '/' + prefix + '_' + sample + '.pbs'
            # test
            with open(newFileName, 'wt') as newFile:
                newFile.writelines(l)
    return 1


def ModifyAndCreate_v2(modelfile, Path_dict, TargetFolder, SampleList, prefix):
    """under each sample folder, create its own PBS files, with name prefix_sample.pbs"""
    with open(modelfile) as file:
        l = file.readlines()
        for i in range(len(l)):
            # Set the job name
            if(l[i].startswith("#PBS -N")):
                l[i] = "#PBS -N " + prefix + '\n'

            elif(l[i].startswith("#PBS -m abe -M")):
                l[i] = "#PBS -m abe -M " + Path_dict.get('Email') + '\n'

            elif(l[i].startswith("wkd=")):
                l[i] = "wkd=" + Path_dict.get('rescue_folder') + '\n'

            elif(l[i].startswith("GATK_Bundle")):
                l[i] = "GATK_Bundle=" + \
                    Path_dict.get("GATK_bundle_path") + '\n'

            elif(l[i].startswith("BWA")):
                l[i] = "BWA=" + Path_dict.get("BWA_path") + '\n'

            elif(l[i].startswith("SAMTOOLS")):
                l[i] = "SAMTOOLS=" + Path_dict.get("samtools_path") + '\n'

            elif(l[i].startswith("PICARD")):
                l[i] = "PICARD=" + Path_dict.get("picard_path") + '\n'

            elif(l[i].startswith("GATK")):
                l[i] = "GATK=" + Path_dict.get("GATK_path") + '\n'

            elif(l[i].startswith("ANNO")):
                l[i] = "ANNO=" + Path_dict.get("ANNO_path") + '\n'

            elif(l[i].startswith("#replace this line by all samples")):
                l[i] = None
                allsampleline = ''
                for sample in SampleList:
                    allsampleline += "    --variant $vcff/" + sample + \
                        '/' + sample + '_RawVariants.g.vcf \\' + '\n'
                l[i] = allsampleline

    newFileName = TargetFolder + '/' + prefix + '.pbs'
    with open(newFileName, 'wt') as newFile:
        newFile.writelines(l)
    return 1
