#######################################################
# Author :  Jainy Thomas
# date   :  December 2019
# email  :  jainythomas1@gmail.com
# Purpose :  extract soft clipped reads from a bam file and blat with query sequencr
#
#####################################################
changelog = '''changelog:
  - v1.0 = Dec 5 2019
    


'''
usage = '''
python3 Fixmate_Transurveyor.py -t Test_SGDP_BAMid.txt -bl /kbod2/WGS_DATA/SGDP_bams_public -o SGDP_outfix 
-pt /home/jainy/software/picard-2.18.23 -pro 5 -st ~/miniconda3/bin 

       '''

import argparse

import subprocess
import os


##################################___LOAD AND CHECK___#####################################
###########################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--table", help="file containing BAMIDs and bamnames")
parser.add_argument("-o", "--outputdir", help="path and name of the output dir")
parser.add_argument("-bl", "--bamlocation", help="location of bam files", type=str)
parser.add_argument("-pro", "--cpu", help="max no. of cpus available")
parser.add_argument("-pt", "--picardtool", help = "path of picard tool")
parser.add_argument("-st", "--samtoolpath", help = "path of samtools")

args = parser.parse_args()

####################################___FUNCTIONS___########################################
###########################################################################################
def check_dir_exists(outdir):
    if os.path.isdir(outdir):
        pass
    else:
        os.makedirs(outdir)
def openfiletoload (filename):
    datadict = {}
    infilobj = open (filename,"r")
    for line in infilobj:#read line by line
        line = line.rstrip()  #
        cols = line.split("\t")
        datadict[cols[0]] = cols[1]

    infilobj.close()
    return(datadict)


def fixmate_picard(picardpath,bampath,listbams,cpus,outdir):
    check_dir_exists(outdir)
    tmpdir = outdir + "/tmp"
    #check_dir_exists(tmpdir)
    child_process = []
    for bamname in listbams:
        bamID =listbams[bamname]
        nbamID = bamID[:-4]
        print("bamID is",bamID, "\n")

        #outBAM = outdir + "/" + nbamID + "_fixed_mate.bam"
        outBAM = outdir + "/" + bamID
        p = subprocess.Popen('java -jar {} FixMateInformation I={} TMP_DIR={} VALIDATION_STRINGENCY=LENIENT O={}'.format(picardpath + "/picard.jar", bampath + "/" + bamID,tmpdir, outBAM), shell=True)
        child_process.append(p)
        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue

    for cp in child_process:
        cp.wait()



def sort_samtools(samtoolspath,listbams,cpus,outdir):
    child_process = []
    for bamname in listbams:
        bamID = listbams[bamname]
        nbamID = bamID[:-4]

        newID = nbamID + "_fixed_mate"
        p = subprocess.Popen(
            '{} sort {} > {}'.format(samtoolspath + "/samtools", outdir + "/" + newID + ".bam", outdir + "/" + newID + "_sorted.bam"), shell=True)
        child_process.append(p)
        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue

    for cp in child_process:
        cp.wait()

def index_samtools(outdir, samtoolspath, listbams, cpus):
    child_process = []
    for bamname in listbams:
        bamID = listbams[bamname]
        nbamID = bamID[:-4]

        #newID = nbamID + "_fixed_mate.bam"

        p = subprocess.Popen(
            '{} index {}'.format(samtoolspath + "/samtools", outdir + "/" + bamID), shell=True)
        child_process.append(p)
        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue

    for cp in child_process:
        cp.wait()



#######################################___MAIN___##########################################
###########################################################################################

bamIDlist = openfiletoload(args.table)

fixmate_picard(args.picardtool, args.bamlocation, bamIDlist, args.cpu, args.outputdir)

#sort_samtools(args.samtoolpath, bamIDlist, args.cpu, args.outputdir)#may be not needed picard tools output in coordsorted manner

index_samtools(args.outputdir, args.samtoolpath, bamIDlist, args.cpu)