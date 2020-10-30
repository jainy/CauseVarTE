#######################################################
# Author :  Jainy Thomas
# date   :  November 2019
# email  :  jainythomas1@gmail.com
# Purpose :  extract soft clipped reads from a bam file and blat with query sequencr
#
#####################################################
changelog = '''changelog:
  - v1.1 = Nov 10 2019
    #added symlink
    #added coverage
    #added MELT-run, Transurveyor,
    #parsed MELT output and Transurveyor output (checked for RM output to see with BP2 overlapps with a TE)
    #Merged both output Run AnnotSV
    #parsed AnnotSV output, ability select based on scores or predicted genotypes
    #Validates insertions by looking for TSDs
    #Rerun Annot
- v1.2 = May 7 2020
    # checked for discordant read mate position intersect with RM data and count the number of occurances of TEs
- v1.3 = June 3 2020
    # removed the discordant read option, as it is time consuming and increases the false positives. reduced the min cut off to 10% reads
- v1.4 = June 3 2020
    # removed the option of identifying TEs from RM data from Transurveyor parse step to reduce time increases the chance for false negatives if you exclude based on that
    # so will have to modify the merge script to exclude the identity for the TE identified
- v1.7 = September 14 2020
    # To do: When I merge, may be keep the MElt genotype as such, (or the complete, still have the Genotype not quality metrics)
- v1.8 = October 7 2020
   # introducing config file and move the arguments from given at the command line
   # Todo to incorporate target site deletions
   # To do, instead of serial numbers provide individual names
   # use parallel for AnnotSV analysis /TSD analysis
    

'''
usage = '''
python3 TEpipeline_script_v1.8.py -conf causevarTE.conf
       '''
import configparser
import argparse
import re
import subprocess
import pysam
import copy
import collections
import os
from decimal import Decimal
import pandas as pd
import shutil
from Bio import AlignIO
from Bio.Align import AlignInfo

##################################___LOAD AND CHECK___#####################################
###########################################################################################

parser = argparse.ArgumentParser()

parser.add_argument("-conf", "--configfile",help = "path of the config file with the parameters")

args = parser.parse_args()

config = configparser.ConfigParser()
config.read(args.configfile)

###accessing toolpaths###########################
musclepath = config.get('Paths','muscle')
AnnSVpath = config.get('Paths','annotSV')
#print(AnnSVpath)
blastpath = config.get('Paths','blast')
meltpath = config.get('Paths','melt')
picardpath = config.get('Paths','picard')
coveragetool = config.get('Paths','coverage')
bwapath = config.get('Paths','bwa')
samtoolspath = config.get('Paths','samtools')
trasuveypath = config.get('Paths','transurveyor')
bedtoolspath = config.get('Paths','bedtools')

###accessing other related info #####################

Refergenome = config.get('refPaths','ref_genome')
bamfiles = config.get('refPaths','bampath')
list_bamid = config.get('refPaths','bamlist')
TE_libr = config.get('refPaths','TElib')
outputdir = config.get('refPaths','outputdir')
hg_ver = config.get('refPaths','hgversion')
RepMaskdata = config.get('refPaths', 'Repeatmaskerdata')
#red
###accessing family info #######################

tot_indi = config.get('famspecific','totaln')
seg_indi = config.get('famspecific','segindn')
nonseg_indi = config.get('famspecific','nonsegind')

###accessing compute info ######################
nupro = config.getint('compute','cpu')
numem = config.getint('compute','mem')

###accessing AnnotSV parsing info ##############
min_score = config.get('AnnotSV','minscore')
colinfo = config.get('AnnotSV','col_info')
colvalue = config.get('AnnotSV','col_valu')

###accessing Processing Conditions##############
consensus = config.getboolean('proconditions','align')
mergelen = config.getint('proconditions','mergedbp')

##########################################################

global TElist
TElist = ["ALU", "SVA", "LINE1"]

# to get the directory containing the script
direname = os.path.dirname(os.path.realpath(__file__))


mini_rd_bpqual = 15 #minimum read base pair quality (average)

####################################___FUNCTIONS___########################################
###########################################################################################
def check_dir_exists(outdir):
    if os.path.isdir(outdir):
        pass
    else:
        os.makedirs(outdir)

def createsimilink(bampath, outputpath, iddict):
    newBamdir = outputpath + "/" + "symBAMdir"
    check_dir_exists(newBamdir)
    # before checking check if already exists or completed

    try:
        with open(outputpath + "/status_check" + "/status_check_file.txt", 'r') as f:
            for line in f:
                line = line.rstrip()
                found = re.search(r"symlink created", line)
                if found:
                    print("SymLINK of bams already created, \n")
    except:
        # create symlink for bam or bai file

        for bamid in iddict:
            bamname = iddict[bamid]
            finalbam = bampath + "/" + bamname
            # newBamdir = newBamdir + "/"
            s = subprocess.run('ln -s {} {} '.format(finalbam, newBamdir), shell=True)

            # print('ln -s {} {} '.format(finalbam, newBamdir))

            # if asmbresultsingle !=0:
            # raise subprocess.CalledProcessError()
            ########in case for .bai does not have .bam not as a prefix

            bamID = bamname[:-4]  # remove .bam from the file name
            print(bamID)
            print(bamname)
            if os.path.exists(bampath + "/" + bamID + ".bai") == True:
                si = subprocess.run('ln -s {} {} '.format(bampath + "/" + bamID + ".bai", newBamdir + "/"), shell=True)
            ########
            elif os.path.exists(bampath + "/" + bamname + ".bai") == True:
                si = subprocess.run('ln -s {} {} '.format(bampath + "/" + bamname + ".bai", newBamdir + "/"), shell=True)
            ########
            # print('ln -s {} {} '.format(bampath + "/" + bamname + ".bai",  newBamdir + "/"))
            # if result !=0:
            # raise subprocess.CalledProcessError()
        with open(outputpath + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
            f_obj.write("symlink created" + "\n")

    return newBamdir

# open file and load the data into a dictionary
def openfiletoload(filename):
    datadict = {}
    infilobj = open(filename, "r")
    for line in infilobj:  # read line by line
        line = line.rstrip()  #
        cols = line.split("\t")
        datadict[cols[0]] = cols[1]

    infilobj.close()
    return (datadict)

def getcoverage(filename):
    datadict = {}
    infilobj = open(filename, "r")
    for line in infilobj:  # read line by line
        if not line.startswith("coverage"):
            line = line.rstrip()  #
            cols = line.split("\t")
            # = cols[12]
            bamname = cols[12].split("/")

            idname = bamname[-1]
            datadict[idname] = cols[0]
    infilobj.close()
    print(datadict)
    return (datadict)

def findcoverage(toolpath, pathdirname, outdir):
    print(toolpath)
    outputfile = outdir + "/bamIDs_coverage.txt"
    print(outputfile)

    if os.path.isfile(outputfile):
        coverage_dict = getcoverage(outputfile)
        return (coverage_dict)
    else:
        asmbresultsingle = subprocess.run(
            '{} covstats {} > {}'.format(toolpath, pathdirname + "/*.bam", outputfile),
            shell=True)

        if asmbresultsingle.returncode != 0:
            raise subprocess.CalledProcessError()
        print('{} covstats {} > {}'.format(toolpath, pathdirname + "/*.bam", outputfile))
        coverage_dict = getcoverage(outputfile)
        return (coverage_dict)





    #return (coverage_dict)
    # print(coverage_dict)
# run MELT(SPLIT)
# preprocess
def run_MELT_SPLIT_preprocess(outdir, melt_path, ref_genome, coverage_dict, bampath, mem, cpus):
    child_process = []  # as part of the parallel through Popen
    for bamname in coverage_dict:
        print(bamname)
        # bam_coverage = coverage_dict[bamname]

        # This is to start the process in parallel using an inbuilt methods
        p = subprocess.Popen(
            'java -Xmx{}g -jar {} Preprocess -bamfile {} -h {}'.format(mem, melt_path + "/MELT.jar",
                                                                       bampath + "/" + bamname,
                                                                       ref_genome), shell=True)
        child_process.append(p)

        # Need to check if the below command work to control the number of jobs launched at a time....

        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue
    for cp in child_process:
        cp.wait()

    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("MELT preprocess completed" + "\n")

    # IndivAnalysis
# Indivanalysis
def run_MELT_SPLIT_indivanal(melt_path, ref_genome, coverage_dict, bampath, mem, cpus, outdir, genomver):
    if genomver.endswith("g38"):
        MEIpath = melt_path + "/me_refs" + "/Hg38"
    elif genomver.endswith("g19"):
        MEIpath = melt_path + "/me_refs" + "/1KGP_Hg19"
    # better to same TE for all Bams so that I can proceed to group analysis once finished
    # calling SVA

    child_process = []  # as part of the parallel through Popen
    for ME in TElist:

        for bamname in coverage_dict:
            print(bamname)
            bam_coverage = coverage_dict[bamname]
            p = subprocess.Popen('java -Xmx{}g -jar {} IndivAnalysis -bamfile {} -h {} -c {} -d {} -w {} -t {}'.format(
                mem, melt_path + "/MELT.jar", bampath + "/" + bamname,
                ref_genome, bam_coverage, "2000000", outdir + "/" + "i" + ME + "dir", MEIpath + "/" + ME + "_MELT.zip"),
                shell=True)
            child_process.append(p)
            if len(child_process) == int(cpus):
                for cp in child_process:
                    cp.wait()
                    child_process = []
            else:
                continue

    for cp in child_process:
        cp.wait()

    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("MELT individual analysis completed" + "\n")
# Group Analysis
def run_MELT_split_groupanal(melt_path, ref_genome, mem, cpus, outdir, genomver):
    if genomver.endswith("g38"):
        MEIpath = melt_path + "/me_refs" + "/Hg38"
        MEIannpath = melt_path + "/add_bed_files" + "/Hg38"
        annfile = "Hg38.genes.bed"
    elif genomver.endswith("g19"):
        MEIpath = melt_path + "/me_refs" + "/1KGP_Hg19"
        MEIannpath = melt_path + "/add_bed_files" + "/1KGP_Hg19"
        annfile = "hg19.genes.bed"

    child_process = []  # as part of the parallel through Popen
    for ME in TElist:

        p = subprocess.Popen('java -Xmx{}g -jar {} GroupAnalysis -h {} -discoverydir {} -w {} -t {} -n {}'.format(
            mem, melt_path + "/MELT.jar",
            ref_genome, outdir + "/" + "i" + ME + "dir", outdir + "/" + "g" + ME + "dir",
                 MEIpath + "/" + ME + "_MELT.zip", MEIannpath + "/" + annfile),
            shell=True)
        child_process.append(p)
        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue
    for cp in child_process:
        cp.wait()

    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("MELT group analysis completed" + "\n")
# Genotype
def run_MELT_split_genotype(melt_path, ref_genome, mem, cpus, outdir, genomver, coverage_dict, bampath):
    if genomver.endswith("g38"):
        MEIpath = melt_path + "/me_refs" + "/Hg38"
        # MEIannpath = melt_path + "/add_bed_files" + "/Hg38"
        # annfile = "Hg38.genes.bed"
    elif genomver.endswith("g19"):
        MEIpath = melt_path + "/me_refs" + "/1KGP_Hg19"
        # MEIannpath = melt_path + "/add_bed_files" + "/1KGP_Hg19"
        # annfile = "hg19.genes.bed"
    print("I am here in genotyping1,\n")
    child_process = []  # as part of the parallel through Popen
    for ME in TElist:
        print("I am here in genotyping2,\n")
        print(ME, "\n")
        for bamname in coverage_dict:
            print(bamname)
            # print('java -Xmx{}g -jar {} Genotype -bamfile {} -h {} -t {} -p {} -w {} '.format(mem, melt_path + "/MELT.jar", bampath + "/" + bamname, ref_genome, MEIpath + "/" + ME + "_MELT.zip", outdir + "/" + "g" + ME + "dir", outdir + "/" + "gt" + ME + "dir"))
            print("I am here in genotyping3,\n")
            # bam_coverage = coverage_dict[bamname]
            p = subprocess.Popen('java -Xmx{}g -jar {} Genotype -bamfile {} -h {} -t {} -p {} -w {} '.format(
                mem, melt_path + "/MELT.jar", bampath + "/" + bamname,
                ref_genome, MEIpath + "/" + ME + "_MELT.zip", outdir + "/" + "g" + ME + "dir",
                     outdir + "/" + "gt" + ME + "dir"),
                shell=True)
            child_process.append(p)
            if len(child_process) == int(cpus):
                for cp in child_process:
                    cp.wait()
                    child_process = []
            else:
                continue

    for cp in child_process:
        cp.wait()

    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("MELT Genotype completed" + "\n")
# MakeVCF
def run_MELT_split_makeVCF(melt_path, ref_genome, mem, cpus, outdir, genomver):
    foutdir = outdir + "/" + "allMEIvcfs"

    check_dir_exists(foutdir)

    if genomver.endswith("g38"):
        MEIpath = melt_path + "/me_refs" + "/Hg38"


    elif genomver.endswith("g19"):
        MEIpath = melt_path + "/me_refs" + "/1KGP_Hg19"

    child_process = []
    for ME in TElist:

        p = subprocess.Popen('java -Xmx{}g -jar {} MakeVCF -ac -genotypingdir {} -h {} -t {} -p {} -w {} -o {}'.format(
            mem, melt_path + "/MELT.jar", outdir + "/" + "gt" + ME + "dir", ref_genome,
                 MEIpath + "/" + ME + "_MELT.zip", outdir + "/" + "g" + ME + "dir", outdir + "/" + "v" + ME + "dir",
                 outdir + "/" + "allMEIvcfs"),
            shell=True)
        child_process.append(p)
        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue
    for cp in child_process:
        cp.wait()

    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("MELT Make VCF completed" + "\n")
# parse MELT

# get breakpoints, get type of TE
# Fixmate neccessary for TS, time consuming step and needs access to hardbam files, recommend to run it separately for the script
def Fixmate_Transur(bamidfile,bampath,picardtoolpath,cpus,outdir,samtoolspath):
    outfixbampth = outdir + "/FixedmateBAMs"
    #tmpdir = outfixbampth + "/tmp"
    check_dir_exists(outfixbampth)
    cmd = "python3 " + direname + "/Fixmate_Transurveyor_v1.0.py " + "-t "+ bamidfile + " -bl "+ bampath + " -pt " + picardtoolpath + " -pro " + cpus + " -o " + outfixbampth + " -st " + samtoolspath
    os.system(cmd)
    #subprocess.call(cmd,shell = True)
    return outfixbampth
     # python3 Fixmate_Transurveyor.py -t Test_BAM_1KGP.txt -bl hg19_1KGP_bams -pt /home/jainy/software/picard-2.18.23 -pro 5 -o hg19_1KGP_bams_mate -st ~/miniconda3/bin > fixmate_status_1.txt &


# Run transurveyor
def Run_transurveyor(tranpath, fixbampath, outdir, refgenome, coverage_dict, cpus, bwapath, samtoolspath):
    print("Running Transurveyor")
    foutdir = outdir + "/" + "Transurveyor_output"
    if os.path.isdir(foutdir):
        shutil.rmtree(foutdir)
        #os.removedirs(foutdir)
        print(foutdir)
    check_dir_exists(foutdir)
    child_process = []
    for bamname in coverage_dict:
        bamID = bamname[:-4]  # remove .bam from the file name
        #nbamID = bamID + "_fixed_mate.bam"
        #print(nbamID)
        bamdir = foutdir + "/" + bamID
        print(bamdir)

        check_dir_exists(bamdir)
        p = subprocess.Popen(
            'python3 {} {} {} {} --threads {} --bwa {} --samtools {} '.format(
                tranpath + "/surveyor_p3_mod.py", fixbampath + "/" + bamname, bamdir, refgenome,
                cpus, bwapath + "/bwa",
                samtoolspath + "/samtools",
            ),
            shell=True)
        #err = p.communicate()[1]#delays the next job to be launched# It also gave the error as none
        #print(err,"HERE I AM ERROR")
        child_process.append(p)
        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue
    for cp in child_process:
        cp.wait()
    if cp.returncode > 0 <= 125:
        print("command failed", cp.returncode)
            #return cp.returncode
        raise subprocess.CalledProcessError()
    elif cp.returncode == 0:
        print("command success", cp.returncode)

            #return cp.returncode

def transurveyor_filter(coverage_dict,tranpath,outdir,cpus):
    print("filtering Transurveyor output")
    foutdir = outdir + "/" + "Transurveyor_output"
    child_process = []
    for bamname in coverage_dict:
        bamID = bamname[:-4]  # remove .bam from the file name
        # nbamID = bamID + "_fixed_mate.bam"
        # print(nbamID)
        bamdir = foutdir + "/" + bamID
        print(bamdir)
        p = subprocess.Popen(
            ' {} {} no-filter > {} '.format(
                tranpath + "/filter", bamdir, bamdir + "/predictions.out"
            ),
            shell=True)
        child_process.append(p)
        if len(child_process) == int(cpus):
            for cp in child_process:
                cp.wait()
                child_process = []
        else:
            continue

    for cp in child_process:
        cp.wait()
    if cp.returncode == 0:

        with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
            f_obj.write("Transurveyor filter completed" + "\n")
    else:
        raise subprocess.CalledProcessError()



    #include check process is completed successfuly or to exit from it

# incorporate ERVcaller identify SVAs in Julie pedigree, which MELT failed but transurveyor identified try with ERVcaller to if it catches,
def print_dict_of_list(TE, dict, outdir):
    with open(outdir + "/allMEIvcfs/allTEconcat.final_comp.bed","a+") as fconct:
        with open(outdir + "/allMEIvcfs" + "/" + TE + ".final_comp.bed", 'w+') as fobj:

            for numbers in sorted(dict.keys()):
                Meltline = dict[numbers]
                fobj.write("\t".join(map(str, Meltline)) + '\n')
                fconct.write("\t".join(map(str, Meltline)) + '\n')

def parse_Melt_VCF_to_bed(outdir, prog):
    print("parsing MELT output")
    for ME in TElist:
        with open(outdir + "/allMEIvcfs" + "/" + ME + ".final_comp.vcf", "r") as f_obj:
            ME_dict = ME + "_dict"
            ME_dict = {}
            num = 1

            for line in f_obj:  # read line by line
                line = line.rstrip()  #

                if line.startswith("##"):
                    continue
                elif line.startswith("#CHROM"):
                    colhead = line.split("\t")
                    totalcol = len(colhead)
                    finalhead = ["#SVchrom", "SV_brk_L", "SV_brk_R", "SV_type", "SV_qual", "TSD", "Strand", "MEinfo"]
                    indivlist = colhead[9:]
                    for indivi in indivlist:
                        finalhead.append(indivi)
                    ME_dict[num] = finalhead
                    #print(finalhead, "title of the bed file\n")
                else:
                    num = num + 1
                    # print(num)
                    bedcol = []

                    cols = line.split("\t")
                    i = 9

                    chrom = cols[0]
                    svbrkpt = cols[1]
                    svbrkpt_L = int(svbrkpt) - 250
                    svbrkpt_R = int(svbrkpt) + 250
                    sv_type = ME
                    sv_qual = prog
                    info = cols[7].split(";")
                    tsd = info[0].split("=")[1]
                    # strand = info[5]
                    MEinfostra = info[5].split("=")[1]
                    strand = MEinfostra[-1]
                    MEinfo = MEinfostra.split(",")[0]

                    bedcol = [chrom, svbrkpt_L, svbrkpt_R, sv_type, sv_qual, tsd, strand, MEinfo]
                    while i < totalcol:
                        indiv = cols[i]
                        bedcol.append(indiv)
                        i = i + 1
                    ME_dict[num] = bedcol
                    # print(bedcol)
            # print(ME_dict)
            print_dict_of_list(ME, ME_dict, outdir)
    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("parsed MELT output" + "\n")
def load_Repeatmasker(pathtofile):
    repeatmaskerdata = {}
    with open(pathtofile, 'r+') as fobj:
        for line in fobj:
            line = line.rstrip()
            reparr = line.split("\t")
            chro = reparr[0]
            start = reparr[1]
            end = reparr[2]
            #Here I am using tuple to use a combine mulitple values to make a unique key for dictionary
            tupkey = (chro,start,end)
            #repeatmaskerdata[chro+":"+start] = []

            # if reparr[4] == "Simple_repeat" or reparr[4] == "Low_complexity":
            #     continue
            if reparr[4] == "Retroposon/SVA" or reparr[4] == "LINE/L1" or reparr[4] == "SINE/Alu":

                #repeatmaskerdata[chro][start] = []#each key has to be unique so this wont work
                #repeatmaskerdata[chro+":"+start] = [chro, start, reparr[2], reparr[3], reparr[4]]
                repeatmaskerdata[tupkey] = [chro, start, reparr[2], reparr[3], reparr[4]]
    #print(repeatmaskerdata)
    return (repeatmaskerdata)
def find_corresponding_TE(chrm,start,end,buffer,rdict):
    TEclass = []
    TEfamil = []
    for chromkey in sorted(rdict.keys()):
        #chrom = chromkey.split(":",0)
        if chromkey[0] == chrm:
            #print(rdict[chromkey]) #eg:'chr21', '9963070', '9963366', 'AluSx1', 'SINE/Alu'

            if int(start) > (int(chromkey[1]) - buffer) and int(start) < int(chromkey[2]) + buffer or \
                    int(end) > (int(chromkey[1]) - buffer) and int(end) < int(chromkey[2]) + buffer:
                # print(rdict[chromkey])
                # print(rdict[chromkey][3])
                # print(rdict[chromkey][4].split("/")[1])
                spr = re.search(r"/",rdict[chromkey][4])
                if spr:
                    transps = rdict[chromkey][4].split("/")[1]
                else:
                    transps = rdict[chromkey][4]
                if transps == "L1":
                    transps = "LINE1"
                elif transps == "Alu":
                    transps = "ALU"
                TEclass.append(transps)
                TEfamil.append(rdict[chromkey][3])
                #print(TEclass, "inside function")
                #print(TEfamil, "inside function")

    return TEclass, TEfamil
def parse_Transurveyor_out(outdir, coverage_dict):
    print("parsing Transurveyor output")
    
    dictallbrkpt = {}
    header = ["#SVchrom", "brkpL", "brkpR", "TEclass", "T", "#discordtreads_Score","inserted_seq", "MEinfo"]
    # Keynumber = range (1,len(coverage_dict))
    # print(Keynumber,"the list created with number of individuals\n")
    nindvi = 1
    for bamname in sorted(coverage_dict.keys()):
        bamID = bamname[:-4]  # remove .bam from the file name
        # nbamID = bamID + "_fixed_mate.bam"
        # print(nbamID)
        header.append(bamID)
        foutdir = outdir + "/" + "Transurveyor_output"
        bamdir = foutdir + "/" + bamID
        print(bamdir)
        #print(nindvi, "the number of individual")
        with open(bamdir + "/predictions.out", "r") as fobj:
            lastline = subprocess.check_output('tail -n 1 {}'.format(bamdir + "/predictions.out"), shell=True)
            lastline = lastline.decode('utf-8')
            # colsp = lastline.split("\t")
            #print(lastline, "lastline")
            IDlist = []
            linenumber = 0
            for line in fobj:  # read line by line
                line = line.rstrip()  #
                cols = line.split(" ")
                #print(line)
                if line == lastline:
                    linenumber = 1
                ID = cols[0]
                #print(IDlist, "Here step1\n")
                if not IDlist:
                    IDlist.append(ID)
                    BP1ar = cols[1].split(":")
                    chr_L = BP1ar[1]
                    brkp_L = BP1ar[2]
                    Key_L = chr_L + ":" + brkp_L
                    BP2ar = cols[2].split(":")
                    BP2chr_L = BP2ar[1]
                    BP2brk_L = BP2ar[2]

                    nDISar = cols[3].split("=")
                    disrds_L = nDISar[1]
                    SCORar = cols[4].split("=")
                    Score_L = round(Decimal(SCORar[1]),1)
                    #print(IDlist, "Here step2\n")
                #when both IDs are equal
                elif IDlist[0] == ID:
                    BP1ar = cols[1].split(":")
                    chr_R = BP1ar[1]
                    brkp_R = BP1ar[2]
                    Key_R = chr_R + ":" + brkp_R
                    BP2ar = cols[2].split(":")
                    BP2chr_R = BP2ar[1]
                    BP2brk_R = BP2ar[2]

                    nDISar = cols[3].split("=")
                    disrds_R = nDISar[1]
                    SCORar = cols[4].split("=")
                    Score_R = round(Decimal(SCORar[1]),1)
                    if (chr_L == chr_R) and abs(int(brkp_L) - int(brkp_R)) < 1000:
                        brkp_F = (int(brkp_L) + int(brkp_R)) / 2
                        brkp_F = int(brkp_F)

                        if int(brkp_L) > int(brkp_R):
                            brkp_Lm = brkp_R
                            brkp_Rm = brkp_L
                        else:
                            brkp_Lm = brkp_L
                            brkp_Rm = brkp_R



                        if Decimal(Score_L) >= Decimal(Score_R):
                            Score_F = Score_L
                        else:
                            Score_F = Score_R

                        if int(disrds_L) >= int(disrds_R):
                            disrds_F = disrds_L
                        else:
                            disrds_F = disrds_R

                        Key_F = chr_R + ":" + str(brkp_F)
                        if BP2chr_L == BP2chr_R:
                            # to change start cordinate always smaller
                            if int(BP2brk_L) < int(BP2brk_R):
                                # TE_class, TE_family = find_corresponding_TE(chr_L, BP2brk_L,BP2brk_R, 100, RMdatadict)
                                #print(TE_class,TE_family,"Here is the TEclass and family")
                                #if not TE_class:
                                TE_class = "NA"
                                TE_family = "NA"
                                # else:
                                #     TE_class = ','.join(TE_class)
                                #     TE_family = ','.join(TE_family)
                                    # TE_class = str(TE_class)[1:-1]# This will keep the ''
                                    #TE_family = str(TE_family)[1:-1]
                                basic_list = [chr_L, brkp_Lm, brkp_Rm, TE_class, "T", disrds_F + "_"+ str(Score_F),
                                          BP2chr_L + ":" + BP2brk_L + "-" + BP2brk_R, TE_family]
                            else:
                                # TE_class, TE_family = find_corresponding_TE(chr_L, BP2brk_R,BP2brk_L, 100, RMdatadict)
                                #if not TE_class:
                                TE_class = "NA"
                                TE_family = "NA"
                                # else:
                                #     TE_class = ','.join(TE_class)
                                #     TE_family = ','.join(TE_family)
                                #print(TE_class,TE_family, "Here is the TEclass and family")
                                basic_list = [chr_L, brkp_Lm, brkp_Rm, TE_class, "T", disrds_F + "_"+ str(Score_F),
                                              BP2chr_L + ":" + BP2brk_R + "-" + BP2brk_L, TE_family]
                        else:
                            # TE_class, TE_family = find_corresponding_TE(chr_L, BP2brk_L, 100, RMdatadict)
                            #if not TE_class:
                            TE_class = "NA"
                            TE_family = "NA"
                            # else:
                            #     TE_class = ','.join(TE_class)
                            #     TE_family = ','.join(TE_family)
                            #print(TE_class,TE_family, "Here is the TEclass and family")
                            basic_list = [chr_L, brkp_Lm, brkp_Rm, TE_class, "T", disrds_F + "_"+ str(Score_F),
                                          BP2chr_L + ":" + BP2brk_L + "_" + BP2chr_R + ":" + BP2brk_R, TE_family]

                        gt = "1/1"
                        nc = "./."
                        for i in range(1, len(coverage_dict) + 1):
                            #print(i)
                            if i == nindvi:
                                basic_list.append(gt)
                            else:
                                basic_list.append(nc)

                        if Key_F in dictallbrkpt.keys():
                            Key_F = Key_F + "_" + str(nindvi)
                            dictallbrkpt[Key_F] = basic_list
                        else:
                            dictallbrkpt[Key_F] = basic_list

                    else:
                        print("breakpoints does not transposition\n")# Todo: write this to a separate file andred were not supported:5000
                    #print(IDlist, "Here step2\n")
                    IDlist = []
                    #print(IDlist, "Here step3\n")
                elif IDlist[0] != ID:



                    brkp_L = int(brkp_L) - 250
                    brkp_R = int(brkp_L) + 250

                    BP2brk_L = int(BP2brk_L) - 50
                    BP2brk_R = int(BP2brk_L) + 50
                    # (TE_class, TE_family) = find_corresponding_TE(chr_L, BP2brk_L, BP2brk_R, 100, RMdatadict)
                    #if not TE_class:
                    TE_class = "NA"
                    TE_family = "NA"
                    # else:
                    #     TE_class = ','.join(TE_class)
                    #     TE_family = ','.join(TE_family)
                    #print(TE_class, TE_family, "Here is the TEclass and family")
                    basic_list2 = [chr_L, brkp_L, brkp_R, TE_class, "T", disrds_L + "_" + str(Score_L),
                                   BP2chr_L + ":" + str(BP2brk_L) + "-" + str(BP2brk_R),TE_family]
                    gt = "1/1"
                    nc = "./."
                    for i in range(1, len(coverage_dict) + 1):
                        #print(i)
                        if i == nindvi:
                            basic_list2.append(gt)
                        else:
                            basic_list2.append(nc)
                    if Key_L in dictallbrkpt.keys():
                        Key_L = Key_L + "_" + str(nindvi)
                        dictallbrkpt[Key_L] = basic_list2
                    else:
                        dictallbrkpt[Key_L] = basic_list2

                    #print(IDlist, "Here step4\n")
                    IDlist = []
                    IDlist.append(ID)
                    BP1ar = cols[1].split(":")
                    chr_L = BP1ar[1]
                    brkp_L = BP1ar[2]
                    key_L = chr_L + ":" + brkp_L
                    BP2ar = cols[2].split(":")
                    BP2chr_L = BP2ar[1]
                    BP2brk_L = BP2ar[2]

                    nDISar = cols[3].split("=")
                    disrds_L = nDISar[1]
                    SCORar = cols[4].split("=")
                    Score_L = round(Decimal(SCORar[1]),1)
                    if linenumber == 1:


                        brkp_L = int(brkp_L) - 250
                        brkp_R = int(brkp_L) + 250

                        BP2brk_L = int(BP2brk_L) - 50
                        BP2brk_R = int(BP2brk_L) + 50
                        # (TE_class, TE_family) = find_corresponding_TE(chr_L, BP2brk_L, BP2brk_R, 100, RMdatadict)
                        #if not TE_class:
                        TE_class = "NA"
                        TE_family = "NA"
                        # else:
                        #     TE_class = ','.join(TE_class)
                        #     TE_family = ','.join(TE_family)
                        #print(TE_class, TE_family, "Here is the TEclass and family")
                        basic_list3 = [chr_L, brkp_L, brkp_R, TE_class, "T", disrds_L + "_"+ str(Score_L),
                                       BP2chr_L + ":" + str(BP2brk_L) + "-" + str(BP2brk_R),TE_family]
                        gt = "1/1"
                        nc = "./."
                        for i in range(1, len(coverage_dict) + 1):
                            # print(i)
                            if i == nindvi:
                                basic_list3.append(gt)
                            else:
                                basic_list3.append(nc)
                        if Key_L in dictallbrkpt.keys():
                            Key_L = Key_L + "_" + str(nindvi)
                            dictallbrkpt[key_L] = basic_list3
                        else:
                            dictallbrkpt[key_L] = basic_list3

                    #print(IDlist, "Here step5\n")
        nindvi += 1
    #print(dictallbrkpt)
    with open(outdir + "/" + "Transurveyor_output/" + "All_TSR_predictions.out", 'w+') as fobj:
        fobj.write("\t".join(header) + '\n')
        for keyele in dictallbrkpt:
            valuele = dictallbrkpt[keyele]
            #print(valuele)
            fobj.write("\t".join(map(str, valuele)) + "\n")
    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("parsed Transurveyor output" + "\n")
def convert_dict_to_dataframe_print(dicta, filename, header,outdir):
    print("converting concatenated, sorted MELT Transurveyor output to pandas")
    print(filename)
    outfilobj = open(filename + ".bef_merge.bed", "w+")
    print("create data from dictionary,\n")
    dfObj = pd.DataFrame.from_dict(dicta, orient='index', columns=["SVLine"])  # creating using dict_Keys as index

    # dfObj = pd.DataFrame(list(dict.items()), index=['a', 'b', 'c', 'd'])
    # print(dfObj)
    newcol = dfObj["SVLine"].str.split("\t", expand=True, )
    #print(header)
    colnames = header.split("\t")
    totcolus = int(len(colnames))
    print(f"total number of columns is {totcolus}")
    #print(colnames)
    newcol.columns = colnames

    #print(newcol)
    #print(newcol.columns)
    # newcol.drop(columns = 0)
    #print(newcol.dtypes)
    newcol2 = newcol.apply(pd.to_numeric, errors='coerce').fillna(newcol)
    #print(newcol2.dtypes)
    newcol2 = newcol2.sort_values(by=[colnames[0], colnames[1]])
    # sortedindex = newcol.sort_index()
    #print(newcol2)
    newcol2[colnames[0]] = newcol2[colnames[0]].astype(str)
    newcol2[colnames[0]] = newcol2[colnames[0]].str.replace('.0', '', regex=False)
    newcol2 = newcol2.sort_values(by=[colnames[0], colnames[1]])
    newcol2.to_csv(outfilobj, sep='\t', index=False)
    outfilobj.close()
    print(" Now Merging the MELT Transurveyor output and print in bed format ")
    mergedfile = filename + ".merge.sorted.bed"
    #print(f"python3 {direname}/Merge_TEcalls_genotypes_v1.3.py -t {filename}.bef_merge.bed -c {totcolus} -l {mergelen} -o {mergedfile}")
    #Todo some close by insertions are merged, needs to modify the script to prevent that or lower the 260 bp to lower
    #or may be I can avoid running ANNOTSV twice run it only once.
    r = subprocess.run(
        'python3 {}/Merge_TEcalls_genotypes_v1.3.py -t {}.bef_merge.bed -c {} -l {} -o {} '.format(
            direname, filename, totcolus,mergelen, mergedfile
        ),
        shell=True)

    if r.returncode == 0:
        with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
            f_obj.write("Merged MELT and Transurveyor output" + "\n")
        return(mergedfile)
    else:
        raise subprocess.CalledProcessError()

    # with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
    #     f_obj.write("Merged MELT and Transurveyor output" + "\n")
    # return(mergedfile)
def load_dict(file,outdir):
    print(" loading the concatenated MELT and Transurveyor output to a dictionary")
    with open(file,"r+") as fo:
        headerarray = []
        allbrkponts = {}
        n = 0
        for line in fo:  # read line by line
            line = line.rstrip()  #
            reparr = line.split("\t")
            if line.startswith("#SVchrom"):
                headerarray.append(line)
                #Todo: need to check if the all columns correctly represents the same individual
            else:
                n = n + 1
                allbrkponts[n] = line
        #print(allbrkponts,"all merged")
    merfedfiletoreturn = convert_dict_to_dataframe_print(allbrkponts, file, headerarray[0],outdir)
    return merfedfiletoreturn
# print(filename, header)
def concat_files(outdir):
    print("concatinating MELT and Transurveyor out")
    concatout = outdir + "/concat_out"
    if os.path.isdir(concatout):
        shutil.rmtree(concatout)
    check_dir_exists(concatout)
    TEfilec = concatout + "/Melt_tsr_concat.txt"
      # concate all output from MELT and transurveyor
    meltcon = outdir + "/allMEIvcfs/allTEconcat.final_comp.bed"#get duplicated list so need to care of it
    TSRout = outdir + "/Transurveyor_output/All_TSR_predictions.out"
    concatfile = subprocess.run('cat {} {} > {} '.format(meltcon, TSRout, TEfilec), shell=True)#
    dict2return = load_dict(TEfilec,outdir)
    return dict2return

def Run_annotSV(AnSVpath,SVinputfile,outdir,genomver,SVminsize):
    # consider that installation requires to the enviroment variable # Todo: automatically set that now the user can set that
    # also exported the path in bash
    print("InsideAnnotsv")
    if genomver == "hg19":
        genomver = "GRCh37"
    elif genomver == "hg38":
        genomver = "GRCh38"
    print(genomver)
    print(SVminsize)
    # if 'ANNOTSV' in os.environ:
    #     print(f"AnnotSV is  set")
    # else:
    #     print(f"AnnotSV is not set")
    os.environ['ANNOTSV'] = AnSVpath
    # if 'ANNOTSV' in os.environ:
    #     print(f"AnnotSV is  set")
    # else:
    #     print(f"AnnotSV is not set")
    AnnotSVout = outdir + "/AnnotSVout"
    check_dir_exists(AnnotSVout)
    inputfilepath = SVinputfile.split("/")[-1]
    #outfilename = SVinputfile[:-4]
    outfilename = inputfilepath[:-4]


    outputffulname = AnnotSVout + "/"+ outfilename + ".annotated.tsv"
    outputallname =  AnnotSVout+"all" + "/"+ outfilename + ".annotated.tsv"
    print(
        'AnnotSV -SVinputFile {} -bedtools {} -outputDir {}  -genomeBuild {} -typeOfAnnotation {} -SVminSize {} -svtBEDcol {}  &> {} '.format(
            SVinputfile, bedtoolspath + "/bedtools", AnnotSVout, genomver, "full", SVminsize, 4,
                         AnnotSVout + "/AnnotSV.log"))

    p = subprocess.run(
        'AnnotSV -SVinputFile {} -bedtools {} -outputDir {}  -genomeBuild {} -typeOfAnnotation {} -SVminSize {} -svtBEDcol {}  &> {} '.format(
            SVinputfile, bedtoolspath +"/bedtools", AnnotSVout,  genomver, "full", SVminsize, 4,  AnnotSVout + "/AnnotSV.log",
        ),
        shell=True)
    q = subprocess.run(
        'AnnotSV -SVinputFile {} -bedtools {} -outputDir {} -genomeBuild {} -SVminSize {} -svtBEDcol {} &> {} '.format(
            SVinputfile, bedtoolspath + "/bedtools", AnnotSVout + "all",  genomver, SVminsize, 4,  AnnotSVout + "/AnnotSV.log",
        ),
        shell=True)
    if p.returncode == 0 and q.returncode == 0:
        with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
            f_obj.write("AnnotSV completed" + "\n")
        return outputffulname, outputallname
    else:
        raise subprocess.CalledProcessError()

    print('AnnotSV -SVinputFile {} -bedtools {} -outputDir {}  -genomeBuild {} -typeOfAnnotation {} -SVminSize {} -svtBEDcol {}  &> {} '.format(
            SVinputfile, bedtoolp +"/bedtools", AnnotSVout,  genomver, "full", SVminsize, 4,  AnnotSVout + "/AnnotSV.log"))

def parse_AnnotSV(annotSVoutfile,totalNumI,scorenu,segrlist,noseglist,outdir,genomver):
    print("inside parse AnnotSV")
    # filenameonly = filename.split("/")[-1]
    #
    # filenamenobed = filenameonly[:-4]  # remove .bed from the file name
    #
    # annotSVoutfile = filenamenobed + ".annotated.tsv"
    # annotSVoutfile
    #parseAnnotSVfile = AnnoSVoutdir +"/"+ filenamenobed +  ".annotated.parse"

    parseAnnotSVfile = annotSVoutfile + "parsed"
    # cmd = "python3 " + direname + "/parseAnnSVTable_forPipeline_v1.0.py -t " + annotSVoutfile + " -n " + totalNumI + " -j 11 -s " + scorenu + " -i " + segrlist + " -a " + noseglist + " -o " + parseAnnotSVfile
    # print(cmd)
    r = subprocess.run('python3 {} -t {} -n {} -j 11 -s {} -i {} -a {} -o {} -vhg {}'.format(direname + "/parseAnnSVTable_forPipeline_v1.0.py",annotSVoutfile,totalNumI,scorenu,segrlist,noseglist,parseAnnotSVfile,genomver),shell= True)
    #os.system(cmd)
    #subprocess.call(cmd,shell = True)
    #returns the input table for validation
    inputfile_discval = parseAnnotSVfile + "_input_table.out"
    inputfile_tsdval = parseAnnotSVfile + "_for_TSDeval.out"
    if r.returncode == 0:
        with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
            f_obj.write("ran and parsed AnnotSV" + "\n")
    else:
        raise subprocess.CalledProcessError()

    return inputfile_discval, inputfile_tsdval
# incorporting TSD identifying script as function
def wholelinetoload (file):
    wholelinedict ={}
    with open (file) as infileobj:
        for line in infileobj:  # read line by line
            line = line.rstrip()  #
            wholelinedict[line]=1
    infileobj.close()
    return(wholelinedict)
def find_avg_list (qual_list):
    return sum(qual_list)/len(qual_list)

def compare_sc_reads(lrdict1,rrdict1,cuoffnu,mincovr):
    print(cuoffnu, "cutoffnumber")
    print(mincovr,"min coverage")
    lrdict1sorbyvalue = OrderedDict(sorted(lrdict1.items(),key=lambda x:x[1],reverse = True))
    rrdict1sorbyvalue = OrderedDict(sorted(rrdict1.items(), key=lambda x: x[1], reverse=True))
    #revsortlrdict1 = sorted(lrdict1.items(),key = lambda item:item[1])# changes to tuple
    #revsortrrdict1 = sorted(rrdict1.items(), key=lambda item: item[1])#changes to tuple
    #tsdcord = "None"
    if len(lrdict1) > 25 or len(rrdict1) > 25:
        print(len(lrdict1))
        print(len(rrdict1))
        return "None"
    else :
        print(len(lrdict1))
        print(len(rrdict1))
        for lk, lv in lrdict1sorbyvalue.items():
            print(lk, lv, "leftclipped")
            for rk, rv in rrdict1sorbyvalue.items():
                print(rk, rv, "rightclipped")
                if cuoffnu <= lv < 200 and cuoffnu <= rv < 200 and 2 < abs(
                        int(rk) - int(lk)) < 26 and int(lv) + int(rv) >= mincovr and rk > lk:  # eventually #9 should be estimated from coverage
                    #last condition would make sure it is a real TSDuplication instead of target site deletion


                    # print("I am here")
                    # numbp = int(rk)-int(lk)
                    # print(numbp)
                    tsdcord = str(lk) + "-" + str(
                        rk)  # add one cordinate to left and substract one right as it 0based cordinate system
                    print("TSDCORd less than 24 is", tsdcord)
                    return tsdcord
                # to catch target site deletions
                # elif cuoffnu <= lv < 200 and cuoffnu <= rv < 200 and 2 < abs(
                #         int(rk) - int(lk)) < 26 and int(lv) + int(rv) >= mincovr and rk < lk:
                #     tsddel = str(rk) + "-" + str(
                #         lk)
                #     print("TSD deletion CORd are", tsddel)
                #     return tsddel # Now there is no distinction target site deletion
        return "None"
                #break

def write_query(sequ, individ, cordid, outdir):
    clippeseq = outdir + "/clippedseqforTEblast" + "/" + individ

    check_dir_exists(clippeseq)
    filetowr = clippeseq + "/" + cordid+ ".combineseq.fasta"
    with open(filetowr, "w+") as quobj:
        quobj.write(">" + cordid + "\n" + str(sequ) + '\n')
    return filetowr

def align_consensus(tsdco,clipseqlist,direction,outdir,indiv,musclepath):
    #print("I am here")
    #print(tsdco)
    muscledir = outdir + "/alignfolder" + "/" +indiv+"/"+tsdco
    inputseqfile = muscledir +"/"+ direction+"clipseq.fasta"
    check_dir_exists(muscledir)

    outfileobj = open(inputseqfile, "w")
    lefcord, rigcord = tsdco.split("-")
    #print(lefcord)
    if direction == "left":
        for elements in clipseqlist:
            cord = elements.split("_")
            tsdcord = cord[0]
            #print(tsdcord)
            cliplen = cord[1]
            if tsdcord == lefcord:
                outfileobj.write(">" + tsdcord + "\n" + cord[2] + '\n')

    elif direction == "right":
        for elements in clipseqlist:
            cord = elements.split("_")
            tsdcord = cord[0]
            cliplen = cord[1]
            if tsdcord == rigcord:
                outfileobj.write(">" + tsdcord + "\n" + cord[2] + '\n')


    outfileobj.close()
    # going align using muscle
    outmusclealn = muscledir +"/"+ direction+"clipseq.align.fasta"
    logfile = muscledir +"/"+ direction+"clipseq.log"
    mresult = subprocess.run(
        '{} -in {} -out {} -log {} -quiet -maxiters {} -diags'.format(musclepath, inputseqfile, outmusclealn, logfile, 3), shell=True)

    # creating consensus using biopython
    alignment = AlignIO.read(outmusclealn,"fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(ambiguous='N',threshold=0.3)
    print(consensus)
    return consensus
def blast_seq (queryseq, dbseq,cpus,Blastpath):

    if os.path.isfile(dbseq + '.nhr'):
        #print("TRUE")
        pass
    else:
        # Keep presets
        BLASTpro = Blastpath + "/" + "makeblastdb"
        result = subprocess.run('{} -in {} -dbtype nucl'.format(BLASTpro,dbseq), shell=True)
        if result.returncode != 0:
            raise subprocess.CalledProcessError()
    BLASTn = Blastpath + "/" + "blastn"
    blastoutput = queryseq + "_blast.out"
    cmd = BLASTn + " -query " + queryseq+ " -db " +dbseq + " -out "+ blastoutput + " -evalue 0.1" + " -num_threads "+ str(cpus) + " -word_size 7" + " -outfmt "+ "\"6 qseqid sseqid pident qlen length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore\""
    #cmd = BLASTn + " -query " + queryseq+ " -db " +dbseq + " -out "+ blastoutput + " -evalue 0.01"

    os.system(cmd)
    return (blastoutput)
    #bresult = subprocess.run('{} -query {} -db {} -out {} -evalue 0.01 -outfmt {}'.format(BLASTn, queryseq, dbseq, blastoutput, "qseqid sseqid pident qlen length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore"), shell=True)
    # if bresult.returncode != 0:
    #     raise subprocess.CalledProcessError()

def longestclippedseq (tsdco,clipseqlist,direction):
    lefcord, rigcord = tsdco.split("-")
    if direction == "left":
        clen =0
        for elements in clipseqlist:
            cord = elements.split("_")
            tsdcord = cord[0]
            #print(tsdcord)
            cliplen = cord[1]
            if tsdcord == lefcord:
                if int(clen) < int(cliplen):
                    clen = cliplen
                    cseq = cord[2]
                else:
                    continue
        return clen, cseq
    elif direction == "right":
        clen = 0
        for elements in clipseqlist:
            cord = elements.split("_")
            tsdcord = cord[0]
            #print(tsdcord)
            cliplen = cord[1]
            if tsdcord == rigcord:
                if int(clen) < int(cliplen):
                    clen = cliplen
                    cseq = cord[2]
                else:
                    continue

        return clen,cseq

def parseblast (blastoutput):
    with open(blastoutput,"r+") as blobj:
        for line in blobj:  # read line by line"qseqid sseqid pident qlen length slen sstrand qstart qend sstart send qcovs mismatch gapopen evalue bitscore"
            line = line.rstrip().split("\t")  #
            predTE = line[1]
            predstrand = line[6]
            #break
            return predTE, predstrand
    return "BND_NoTEinfo","."
def print_finalout(outdir,filename,dict):
    check_dir_exists(outdir)
    with open(outdir +"/"+filename,"w+") as fobj:
        for baid in dict:
            #print(baid)
            # print(dict[baid])
            for uid in dict[baid]:
                #print(uid)
                #print(dict[baid][uid])
                if not dict[baid][uid]["TSDcor"] == "None":
                    fobj.write(baid+"\t"+uid+"\t"+dict[baid][uid]["TSDcor"]+"\t"+dict[baid][uid]["TEpred"]+"\t"+dict[baid][uid]["strand"]+"\t"+dict[baid][uid]['TE5p']+"\t"+dict[baid][uid]['TE3p']+"\n")
                else:
                    fobj.write(baid + "\t" + uid + "\t" + dict[baid][uid]["TSDcor"]+"\n")

def print_finalout_hori(outdir,filename,dict,headerf,totalnumb):
    check_dir_exists(outdir)
    outfile_annotSV = outdir +"/"+filename
    indi =[]
    infocol = list(range(4,9))
    with open(outfile_annotSV,"w+") as fobj:

        fobj.write("\t".join(headerf)+"\n")
        gt = "1/1"
        nc = "./."
        totcol = 8+int(totalnumb) # from 9th column,  genotyping info starts

        #print(indi)
        print(infocol)

        for tsdid in dict:
            print(tsdid)
            wholeline = []
            if len(dict[tsdid]) <= totcol - 4:
                print(len(dict[tsdid]), "length of dict")
                print(totcol, "total column number")
                for k in infocol:
                    if k in dict[tsdid].keys():
                        pass
                    else:
                        dict[tsdid][k] = "."


            indi = list(range(9, int(totalnumb) + 9))
            print(indi)
            if len(dict[tsdid]) < totcol:
                print(tsdid,len(dict[tsdid]))

                for i in indi:
                    if i in dict[tsdid].keys():
                        pass
                    else:
                        dict[tsdid][i] = nc



            for uid in sorted(dict[tsdid].items(),key=lambda x: int(x[0])):
                print(uid)
                wholeline.append(uid[1])


            #print(wholeline)
            fobj.write("\t".join(wholeline)+"\n")
    return outfile_annotSV
def loca_te_discmate(dmatedict,rdict):
    buffer = 0
    leng_dm = len(dmatedict)
    #print(dmatedict)
    print ("dismate dict length",leng_dm)
    supportedTEfamily = {}
    dmdict1sorbyvalue = OrderedDict(sorted(dmatedict.items(), key=lambda x: x[1], reverse=True))
    for loci in dmdict1sorbyvalue:
        #print(loci,dmdict1sorbyvalue[loci])
        chrloci = loci.split("|")
        chrm,start,end = chrloci[0].split(";")


        TEclassy, TEfamily = find_corresponding_TE(chrm, start, end, buffer, rdict)

        if len(TEclassy) ==0 and len(TEfamily) == 0:
            continue
        else:
            #print(TEclassy[0], TEfamily[0], "from find corresponding TE")
            if TEclassy[0] in supportedTEfamily:
                supportedTEfamily[TEclassy[0]] += dmdict1sorbyvalue[loci]
            else:
                supportedTEfamily[TEclassy[0]] = dmdict1sorbyvalue[loci]


    if not len(supportedTEfamily) ==0:
        supportedTEfamilysorted = OrderedDict(sorted(supportedTEfamily.items(), key=lambda x: x[1], reverse=True))
        #print(supportedTEfamilysorted,"TEfamily for each loci")
        for TEelme in supportedTEfamilysorted:
            print(TEelme,supportedTEfamilysorted[TEelme])
            return TEelme,supportedTEfamilysorted[TEelme]
    else:
        return "NoTE", 0


def call_TSD_to_identifybrkpt(anSVout1,coverage_dict, bampath,musclepath,TElibpath,outdir,cpus,blastpath,totalNumI):
    print("I am here within tsd")
    coordstoanalyse = wholelinetoload(anSVout1)
    RMdata_dict = load_Repeatmasker(RepMaskdata)
    #print(coordstoanalyse)
    TSD_cord={}
    gtdata = {}#unique_id Chr1 TSD_start TSD_end strand TE5p_seq TE3q_seq Todo:No.of.discordantreadswithmate_mappedto_otherTEs gtIndi1 gtind2 gtind3 gtind4
    headerforcombgttb = ["#Chr", "TSD_start","TSD_end", "TE","strand","TE3q_seq","TE5p_seq", "TEfamily"]
    bindex = 9#hard core value, this is the column number where individual genotypes starts check if that changes
    TSD_status ={}
    DM_info_TE = {}  # to store the final info from DM read analysis
    for bamfilename in sorted(coverage_dict.keys()):
        coverage = coverage_dict[bamfilename]
        cutoffnumber = round(.08 * float(coverage))
        mincover = round(.10 * float(coverage))
        maxcover = round(2 * float(coverage))
        bam = bamfilename[:-4]
        print(bam)
        print(bindex,"bamindex")

        #All_discmate[bindex]= {}# storing output from the intersect data
        Discmateread ={}# loading the position of discordant mate reads for each unique id in a bam and will sort it and print that to file, Then intersect
        #with RM data , parse the output to find the likely TE inserted and store than in a dictionoary
        headerforcombgttb.append(bam)
        bamfiletopen = bampath + "/" + bamfilename
        #print(bamfiletopen)
        samfile = pysam.AlignmentFile(bamfiletopen, "rb")
        TSD_cord[bam] = {}
        for cords in coordstoanalyse:
            print(cords)
            if cords.startswith("#"):
                continue
            else:
                cols = cords.split( )  # split on space
                # cols = cords.split()
                print(cols[0])
                chr = cols[0]
                # if genomver == "hg19":
                #     chr = int(chr)
                start = cols[1]
                end = cols[2]
                # TE = cols[3]
                svlength = int(end) - int(start)
                uniqueid = chr + ":" + start + "-" + end
                TSD_status[uniqueid] = {}

                if svlength < 250:
                    meanbkp = (int(end) + int(start))/2
                    start = meanbkp -250
                    end = meanbkp +250
                #VCF_gtdata[uniqueid]=[]
                TSD_cord[bam][uniqueid] = {}
                RighSoftclipcord = {}
                LeftSoftclipcord = {}
                RighSoftclipseq = []
                LeftSoftclipseq = []
                Discmateread = {}
                print("Now analysing", uniqueid)
                n = 1
                for read in samfile.fetch(chr, int(start), int(end)):
                    # mapstatus = 3
                    #print("all reads are ",read)

                    if read.is_qcfail == True or read.is_duplicate == True:
                        # print("quality failed reads are", read)
                        continue
                    else:
                        refchr = read.reference_name
                        refstart = read.reference_start
                        refend = read.reference_end
                        readnameonly = read.query_name
                        readseq = read.query_sequence
                        readmapq = read.mapping_quality

                        #print(readseq)
                        # readqaln = read.query_alignment_qualities
                        if read.is_read1 == True:
                            freadname = readnameonly + "/1"
                            matereadname = readnameonly + "/2"
                        elif read.is_read2 == True:
                            freadname = readnameonly + "/2"
                            matereadname = readnameonly + "/1"
                        else:
                            freadname = readnameonly
                            matereadname = readnameonly



                        if read.is_proper_pair == True:
                            typofMap = "proper"
                            flen = read.template_length
                        else:
                            typofMap = "discordant"
                            #print(typofMap)
                            flen = read.template_length
                            nerefchr = read.next_reference_name
                            nerefstart = read.next_reference_start
                            nerefend = int(nerefstart) + int(read.query_alignment_length)
                            #nerefrd = read.next_reference_id
                            #nerefrd, read.reference_id
                            #print(nerefchr,nerefstart,nerefend)
                            dmateid = nerefchr + ";" + str(nerefstart) + ";" + str(nerefend) + "|"+uniqueid

                            if dmateid in Discmateread:
                                Discmateread[dmateid] += 1
                            else:
                                Discmateread[dmateid] = 1

                            n += 1


                        cstring = read.cigarstring
                        #print(cstring)# todo: can do 1 set of filtering here
                        # readprop = pysam.AlignedSegment()
                        tup_list = read.cigartuples
                        #print("cigartuples", read.cigartuples, freadname)
                        if read.cigartuples:
                            lastindex = (len(tup_list) - 1)

                            if lastindex == 0:
                                #print(lastindex, "The last index")
                                continue
                            else:
                                indexofsof = [i for i, y in enumerate(tup_list) if
                                              y[0] == 4 and y[1] > 7]  # len of the softclipped reads more than 20
                                #print(indexofsof)
                                if indexofsof == [0]:
                                    #print("clipped reads",tup_list[0][1])
                                    readstart = read.reference_start
                                    #print("refstart", readstart)
                                    numof_LeftSoftclp = int(tup_list[0][1])
                                    #print("left clipped reads",numof_LeftSoftclp)
                                    lsqua = read.query_qualities[1:numof_LeftSoftclp]
                                    #print(lsqua)
                                    avglqsc = find_avg_list(lsqua)
                                    #print(round(avglqsc))
                                    if avglqsc > mini_rd_bpqual:
                                        lscseq = read.query_sequence[1:numof_LeftSoftclp]
                                        #print(lscseq)
                                        lscseqfas = str(readstart) + "_" + str(numof_LeftSoftclp) + "_" + lscseq
                                        #print(lscseqfas)
                                        LeftSoftclipseq.append(lscseqfas)
                                        if readstart in LeftSoftclipcord:
                                            LeftSoftclipcord[readstart] += 1
                                        else:
                                            LeftSoftclipcord[readstart] = 1

                                    else:
                                        continue


                                if indexofsof == [lastindex]:
                                    numof_RightSoftclp = int(tup_list[lastindex][1])
                                    #print("clippedreads", tup_list[lastindex][1])
                                    readend = read.reference_end
                                    #print("refend", readend)
                                    rsqua = read.query_qualities[-numof_RightSoftclp:]
                                    #print(rsqua)
                                    avgrqsc = find_avg_list(rsqua)

                                    if avgrqsc > mini_rd_bpqual:
                                        # print(round(avgrqsc))
                                        rscseq = read.query_sequence[-numof_RightSoftclp:]
                                        # print(rscseq)
                                        rscseqfas = str(readend) + "_" + str(numof_RightSoftclp) + "_" + rscseq
                                        RighSoftclipseq.append(rscseqfas)

                                        if readend in RighSoftclipcord:
                                            RighSoftclipcord[readend] += 1

                                        else:
                                            RighSoftclipcord[readend] = 1
                                    else:
                                        continue

                        else:
                            continue
                            # if read.is_unmapped == True:
                            #     unmappedrds[bam][freadname] = 1
                            #     unmappedreads[bam][freadname] = readseq
                            #     # print(read)
                            # else:
                            # print("read mapped", read)
                #print(LeftSoftclipcord,"Leftsoftcordinates")
                TSDcord = compare_sc_reads(LeftSoftclipcord, RighSoftclipcord,cutoffnumber,mincover)
                TSD_cord[bam][uniqueid]["TSDcor"] = TSDcord

                # TSD_id = chr + ":" + TSDcord  # adding chromosome
                # TSDstart, TSDend = TSDcord.split("-")
                #   1 = chr
                #   2 = TSDstart
                #   3 = TSDend
                #   4 = TEpred
                #   5 = strand

                #   6 = TE5p

                #   7 = TE3p
                #   8 = TEfamily
                #   9 = genotype (to however many) # Now only 1/1 or ./.
                #
                #

                if not TSDcord == "None":
                    print(bam, TSDcord, "TSD identified")
                    TSD_id = chr + ":" + TSDcord  # adding chromosome
                    TSDstart, TSDend = TSDcord.split("-")
                    if uniqueid in TSD_status:
                        if "infoTSD" in TSD_status[uniqueid].keys():
                            existing_TSDident = TSD_status[uniqueid]["infoTSD"]
                            if existing_TSDident == "None":
                                TSD_status[uniqueid]["infoTSD"] = TSD_id
                                #pass
                        else:
                            TSD_status[uniqueid]["infoTSD"] = TSD_id
                            #print(TSD_status[uniqueid])# Need to compare process later#choose a majority one probability

                    else:
                        TSD_status[uniqueid]["infoTSD"] = TSD_id


                    if TSD_id in gtdata:
                        pass
                    else:
                        gtdata[TSD_id] = {}
                        gtdata[TSD_id][1] = chr
                        gtdata[TSD_id][2] = TSDstart
                        gtdata[TSD_id][3] = TSDend
                        #gtdata[TSD_id]["genotypes"] = []
                    if consensus == True:

                        leftconsensus = align_consensus(TSDcord, LeftSoftclipseq, "left", outputdir, bam,
                                                        musclepath)
                        rightconsensus = align_consensus(TSDcord, RighSoftclipseq, "right", outputdir, bam,
                                                         musclepath)
                        longcombineseq = rightconsensus + leftconsensus
                        print(longcombineseq)
                        TSD_cord[bam][uniqueid]["TE3p"] = str(leftconsensus)
                        TSD_cord[bam][uniqueid]["TE5p"] = str(rightconsensus)
                        # longcombineseq = rightconsensus + leftconsensus
                        lenleftcons = len(leftconsensus)
                        lenrigtcons = len(rightconsensus)

                        if TSD_id in gtdata:
                            #if "TE3p" in gtdata[TSD_id].keys():
                            if 7 in gtdata[TSD_id].keys():
                                existingleftconsensus = gtdata[TSD_id][7]
                                if len(existingleftconsensus) > lenleftcons:
                                    pass
                                else:
                                    gtdata[TSD_id][7] = str(leftconsensus)
                            else:
                                gtdata[TSD_id][7] = str(leftconsensus)

                            #if "TE5p" in gtdata[TSD_id].keys():
                            if 6 in gtdata[TSD_id].keys():
                                existingrightconsensus = gtdata[TSD_id][6]
                                if len(existingrightconsensus) > lenrigtcons:
                                    pass
                                else:
                                    gtdata[TSD_id][6] = str(rightconsensus)


                            else:
                                gtdata[TSD_id][6] = str(rightconsensus)



                        else:
                            gtdata[TSD_id][7] = str(leftconsensus)
                            gtdata[TSD_id][6] = str(rightconsensus)

                        

                    else:
                        maxl, maxlseq = longestclippedseq(TSDcord, LeftSoftclipseq, "left")
                        print(maxl, maxlseq)

                        maxr, maxrseq = longestclippedseq(TSDcord, RighSoftclipseq, "right")
                        print(maxr, maxrseq)
                        TSD_cord[bam][uniqueid]["TE3p"] = maxlseq
                        TSD_cord[bam][uniqueid]["TE5p"] = maxrseq
                        longcombineseq = maxrseq + maxlseq

                    queryfile = write_query(longcombineseq, bam, uniqueid, outdir)
                    # Blast to TE library using combined seq to identify TE inserted and strand of insertion
                    # TEssequne, consensus AluY, LINE-1 and SVA_DEF, HERVK, (HERVW)
                    blastout = blast_seq(queryfile, TElibpath, cpus, blastpath)
                    TEtype, strand = parseblast(blastout)
                    #print(TEtype,strand)
                    TE,TEfamily = TEtype.split("_")
                    #print(TE,TEfamily)
                    TSD_cord[bam][uniqueid]["TEpred"] = TEtype
                    TSD_cord[bam][uniqueid]["strand"] = strand

                    #if "TEpred" in gtdata[TSD_id].keys():
                    if 4 in gtdata[TSD_id].keys():
                        
                        currenTEtype = gtdata[TSD_id][4]
                        if currenTEtype == TE:
                            #print("no change")
                            pass
                        else:
                            #possiTE = currenTEtype + ";" + TEtype
                            possiTE = "MEI"
                            gtdata[TSD_id][4] = possiTE
                            print(gtdata[TSD_id][4])  # may be combine both?

                    else:
                        gtdata[TSD_id][4] = TE
                        
                    if 8 in gtdata[TSD_id].keys():
                        
                        currenTEfamily = gtdata[TSD_id][8]
                        if currenTEfamily == TEfamily:
                            #print("no change")
                            pass
                        else:
                            possiTEfamily = currenTEfamily + ";" + TEfamily
                            #possiTE = "MEI"
                            gtdata[TSD_id][8] = possiTEfamily
                            print(gtdata[TSD_id][8])  # may be combine both?

                    else:
                        gtdata[TSD_id][8] = TEfamily

                    if 5 in gtdata[TSD_id].keys():
                            # existingleftconsensus = gtdata[TSD_id]["TE3p"]
                        currenTEstrand = gtdata[TSD_id][5]
                        if currenTEstrand == strand:
                            #print("no change")
                            pass
                        else:
                            possistrand = currenTEstrand + ";" + strand
                            gtdata[TSD_id][5] = possistrand
                            print(gtdata[TSD_id][5])  # may be combine both?

                    else:

                        gtdata[TSD_id][5] = strand
                    #gtdata[TSD_id]["genotypes"].append("1/1")
                    if bindex in gtdata[TSD_id].keys():
                        currentgt = gtdata[TSD_id][bindex]
                        if currentgt == "./.":
                            gtdata[TSD_id][bindex] = "1/1"

                            print(currentgt, "current genotype")
                    else:
                        gtdata[TSD_id][bindex] = "1/1"
                else:
                    #All_discmate[uniqueid] = "None"
                    if uniqueid in TSD_status:
                        if "infoTSD" in TSD_status[uniqueid].keys():
                            pass
                        else:
                            TSD_status[uniqueid]["infoTSD"] = "None"
                            #print(TSD_status[uniqueid])# Need to compare process later#choose a majority one probability



        bindex += 1
        samfile.close()
        #intersectBed_RMdata(Discmateread,RMdata_dict)

        # print(TSD_cord)
        # print(RighSoftclipseq)
        # print(LeftSoftclipseq)

    #print(TSD_cord)
    #print(gtdata)
    #print(headerforcombgttb)
    print_finalout(outputdir, "refinedtsd.txt", TSD_cord)
    #print(TSD_status,"TSDrecord-uniqueid dict")
    fileforASVrd2 = print_finalout_hori(outputdir, "filetorunAnSV_2.bed", gtdata, headerforcombgttb,totalNumI)
    print(fileforASVrd2, "filename for second AnnotSV")
    with open(outdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
        f_obj.write("TSD call completed" + "\n")

    return fileforASVrd2 



#######################################___MAIN___##########################################
###########################################################################################
# open file and load the data into a dictionary
Bamid_dict = openfiletoload(list_bamid)
print(Bamid_dict)
# create a directory
check_dir_exists(outputdir)
status_checkdir = outputdir + "/status_check"
check_dir_exists(status_checkdir)  # create a dir for status check
SoftBamDir = createsimilink(bamfiles, outputdir, Bamid_dict)
cov_dict = []

if coveragetool.endswith("coverage.txt"):
    cov_dict = getcoverage(coveragetool)  # coverage is stored under whole bamname#
    # print(cov_dict)

else:
    # find coverage for each sample in teh path of the bamfile folder using gostats

    cov_dict = findcoverage(coveragetool, SoftBamDir, outputdir)

print("Coverage dict", cov_dict)
# Run MELT preprocess

try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()

            found = re.search(r"MELT preprocess completed", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
    if count == 1 :
        print("MELT preprocess of bams already done, \n")
    else:
        print("I am here in step1,\n")
        run_MELT_SPLIT_preprocess(outputdir, meltpath, Refergenome, cov_dict, SoftBamDir, numem,
                                  nupro)
except:
    print(
        "the status file is not found. It should be created with the creation of the symlink. However proceeding to MELT-processing...")
    exit()

    # run_MELT_SPLIT_preprocess(outputdir, meltpath, Refergenome, cov_dict, SoftBamDir, numem, nupro)
# Run MELT Indivanalysis
#run_MELT_SPLIT_preprocess(outputdir, meltpath, Refergenome, cov_dict, SoftBamDir, numem,
                                  #nupro)
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"MELT individual analysis completed", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:
            print("MELT individual analysis of bams already done, \n")
        else:
            print("I am here in step2,\n")
            run_MELT_SPLIT_indivanal(meltpath, Refergenome, cov_dict, SoftBamDir, numem, nupro,
                                     outputdir, hg_ver)
except:
    print(
        "the status file is not found or failed to run individual analysis....")
    exit()
    # run_MELT_SPLIT_indivanal(meltpath, Refergenome, cov_dict, SoftBamDir, numem, nupro,
    # outputdir, hg_ver)

##Run MELT Groupanalysis
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"MELT group analysis completed", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("MELT group analysis of bams already done, \n")
        else:
            print("I am here in step3,\n")
            run_MELT_split_groupanal(meltpath, Refergenome, numem, nupro, outputdir,
                                     hg_ver)
except:
    print(
        "the status file is not found or group analysis could not be performed,\n  ")
    exit()

    # run_MELT_split_groupanal(meltpath, Refergenome, numem, nupro, outputdir, hg_ver)

##Run MELT Genotype
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"MELT Genotype completed", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("MELT genotyping already done, \n")
        else:
            print("I am here in step4,\n")
            run_MELT_split_genotype(meltpath, Refergenome, numem, nupro, outputdir, hg_ver,
                                    cov_dict, SoftBamDir)
except:
    print(
        "the status file is not found or genotyping step cannot be performed,\n ")
    exit()

    # run_MELT_split_genotype(meltpath, Refergenome, numem, nupro, outputdir, hg_ver)

##Run MELT makeVCF
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"MELT Make VCF completed", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("MELT final VCFs already made, \n")
        else:
            print("I am here in step5,\n")
            run_MELT_split_makeVCF(meltpath, Refergenome, numem, nupro, outputdir, hg_ver)
except:
    print("the status file is not found or MakeVCF step cannot be performed,\n ")
    exit()

# Run Transurveyor
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"Transurveyor filter completed", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("Transurveyor filter completed already, \n")
        else:
            try:
                print("I am here in step6")  # workspace needs to be emptied then only Transurveyor runs

                Run_transurveyor(trasuveypath, SoftBamDir, outputdir, Refergenome, cov_dict, nupro, bwapath,
                             samtoolspath)
                transurveyor_filter(cov_dict, trasuveypath, outputdir, nupro)
            except:
                print("I am here in step 6.1")
                fix_bampath = Fixmate_Transur(list_bamid, SoftBamDir, picardpath, nupro, outputdir, samtoolspath)
                fix_bampath = outputdir + "/FixedmateBAMs"
                Run_transurveyor(trasuveypath, fix_bampath, outputdir, Refergenome, cov_dict, nupro, bwapath,
                             samtoolspath)
                transurveyor_filter(cov_dict, trasuveypath, outputdir, nupro)


except:
    print(
        "the status file is not found or Transurveyor failed cannot be performed,\n ")
    exit()
####Parsing MELT outpu######
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"parsed MELT output", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("parsed MELT output already, \n")
        else:
            print("I am here in step7,\n")
            parse_Melt_VCF_to_bed(outputdir, "M")
except:
    print(
        "the status file is not found or parse MELT step cannot be performed,\n ")
###parsing Transurveyor output#####

try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"parsed Transurveyor output", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("parsed Transurveyor output already, \n")
        else:
            print("I am here in step8,\n")
            #RMdata_dict = load_Repeatmasker(RepMaskdata)# loaded repeatmasker to annotate the repeats that are originated from TEs
            #print("loaded RM data")
            parse_Transurveyor_out(outputdir, cov_dict)
            #
except:
    print(
        "the status file is not found or parse Transurveyor step cannot be performed,\n ")

# concatenate Melt and Transureveryor,load into dictionary, covert to pandas df, Merge them, and sort  them and print the sorted file
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"Merged MELT and Transurveyor output", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("Merged MELT and Transurveyor output already, \n")
        else:
            print("I am here in step9,\n")
            mergedsortedfile = concat_files(
                outputdir)  # files get bigger all the time because concatenation so need to remove if repeated more than once
            print(mergedsortedfile)
            mergedfilecopy= outputdir + "/mergedfileforAnnotSV.bed"
            shutil.copy(mergedsortedfile, outputdir + "/mergedfileforAnnotSV.bed" )
            #os.system('cp {} {}'.format(mergedsortedfile, mergedfilecopy))

except:
    print(
        "the status file is not found or  Merge step cannot be performed,\n ")
    exit()




#Run ANNOTSV for all merged breakpoints before doing TSD analysis######





try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"ran and parsed AnnotSV", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("AnnotSV completed and parsed already, \n")
        else:

            print("I am here in step10")
            try:
                print("I am here in step10.1")
                mergedsortedfile
                annotfullfile,annotallfile = Run_annotSV(AnnSVpath, mergedsortedfile, outputdir,hg_ver,50)
                print(f"AnnotSV Step 10.1 done {annotfullfile} generated")
                fulinputfileforcord_indv, fulinputfileforonlycord = parse_AnnotSV(annotfullfile, tot_indi,
                                                                                  min_score,
                                                                                  seg_indi,
                                                                                  nonseg_indi, outputdir,hg_ver)
                copyinputtable1 = outputdir + "/fullinputtable_ind_with_brkpt.bed"
                copyinputtable2 = outputdir + "/fullinputtable_only_brkpt.bed"
                shutil.copy(fulinputfileforcord_indv, copyinputtable1)
                shutil.copy(fulinputfileforonlycord, copyinputtable2)
                print(f"parse AnnotSV out Step 10.1 done {fulinputfileforcord_indv} generated")
                #os.system('cp {} {}'.format(inputfileforval, copyinputtable))
            except:
                print("I am here in step10.2")
                mergedsortedfile = outputdir + "/mergedfileforAnnotSV.bed"
                annotfullfile,annotallfile = Run_annotSV(AnnSVpath, mergedsortedfile, outputdir,hg_ver,50)
                print(f"AnnotSV Step 10.2 done {annotfullfile} generated")
                fulinputfileforcord_indv, fulinputfileforonlycord = parse_AnnotSV(annotfullfile, tot_indi,
                                                                                  min_score,
                                                                                  seg_indi,
                                                                                  nonseg_indi, outputdir,hg_ver)
                copyinputtable1 = outputdir + "/fullinputtable_ind_with_brkpt.bed"
                copyinputtable2 = outputdir + "/fullinputtable_only_brkpt.bed"
                shutil.copy(fulinputfileforcord_indv, copyinputtable1)
                shutil.copy(fulinputfileforonlycord, copyinputtable2)
                print(f"parse AnnotSV out Step 10.2 done {fulinputfileforcord_indv} generated")
                #os.system('cp {} {}'.format(inputfileforval, copyinputtable))

except:
    print(
        "the status file is not found or AnnotSV step cannot be performed,\n ")
    exit()


#######Call TSD#####

try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"TSD call completed", line)
            if found:
                count = 1
            else:
                continue
            # print(found)
        if count == 1:

            print("TSD call completed, \n")
        else:
            print("I am here in step11")

            try:
                fulinputfileforonlycord
                print(fulinputfileforonlycord)
                print("I am here in step11a")
                filefor2ndAnnotSV = call_TSD_to_identifybrkpt(fulinputfileforonlycord, cov_dict, SoftBamDir, musclepath, TElib,
                                                                outputdir, nupro, blastpath,tot_indi)
                # annotfulltsd, annotalltsdfile = Run_annotSV(AnnSVpath, filefor2ndAnnotSV, outputdir,
                #                                             hg_ver, "1")
            except:
                print("I am here in step11b")
                fulinputfileforonlycord = outputdir + "/fullinputtable_only_brkpt.bed"

                filefor2ndAnnotSV = call_TSD_to_identifybrkpt(fulinputfileforonlycord, cov_dict, SoftBamDir, musclepath, TElib,
                                                      outputdir, nupro, blastpath,tot_indi)

                # annotfulltsd, annotalltsdfile = Run_annotSV(AnnSVpath, filefor2ndAnnotSV, outputdir,
                #                                             hg_ver, "1")



except:
    print(
        "the status file is not found or TSD finding cannot be performed,\n ")
    exit()

##Perform AnnotSV###

try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"AnnotSV completed on TSDcandidates", line)
            if found:
                count = 1
            else:
                continue
            # print(found)
        if count == 1:

            print("Ran AnnotSV on TSD candidates, \n")
        else:
            print("I am here in step12")

            try:
                filefor2ndAnnotSV
                print(filefor2ndAnnotSV)
                print("I am here in step12a")

                annotfulltsd, annotalltsdfile = Run_annotSV(AnnSVpath, filefor2ndAnnotSV, outputdir,
                                                            hg_ver, "1")
                if (annotfulltsd and annotalltsdfile):
                    with open(outputdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
                        f_obj.write("AnnotSV completed on TSDcandidates" + "\n")

            except:
                print("I am here in step12b")
                filefor2ndAnnotSV = outputdir + "/filetorunAnSV_2.bed"

                annotfulltsd, annotalltsdfile = Run_annotSV(AnnSVpath, filefor2ndAnnotSV, outputdir,
                                                            hg_ver, "1")
                if (annotfulltsd and annotalltsdfile):
                    with open(outputdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
                        f_obj.write("AnnotSV completed on TSDcandidates" + "\n")



except:
    print(
        "the status file is not found or ANNOTSV on TSDcandidates cannot be performed,\n ")
    exit()
####Parse ANNOTSV#####
try:
    with open(outputdir + "/status_check" + "/status_check_file.txt", 'r') as f:
        count = 0
        for line in f:
            line = line.rstrip()
            found = re.search(r"parsed 2nd_AnnotSV output", line)
            if found:
                count = 1
                break
            else:
                continue
            # print(found)
        if count == 1:

            print("parsed final AnnotSV output, \n")
        else:
            print("I am here in step13")

            try:
                annotalltsdfile
                print(annotalltsdfile)
                print("I am here in step13a")
                r = subprocess.run(
                    'python3 {} -t {} -n {} -j 11 -s {} -i {} -a {} -cinfo {} -cval {}'.format(
                        direname + "/parseAnnSVTable_forTSDout_v1.1.py", annotalltsdfile, tot_indi, min_score,
                        seg_indi,
                        nonseg_indi, colinfo, colvalue), shell=True)
                # os.system(cmd)

                if r.returncode == 0:
                    print("step13a success")
                    with open(outputdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
                        f_obj.write("parsed 2nd_AnnotSV output" + "\n")
                else:
                    raise subprocess.CalledProcessError()


            except:
                print("I am here in step13b")
                annotalltsdfile = outputdir + "/AnnotSVoutall/filetorunAnSV_2.annotated.tsv"

                r = subprocess.run(
                    'python3 {} -t {} -n {} -j 11 -s {} -i {} -a {} -cinfo {} -cval {}'.format(
                        direname + "/parseAnnSVTable_forTSDout_v1.1.py", annotalltsdfile, tot_indi, min_score,
                        seg_indi,
                        nonseg_indi, colinfo, colvalue), shell=True)
                if r.returncode == 0:
                    print("step13b success")
                    with open(outputdir + "/status_check" + "/status_check_file.txt", "a+") as f_obj:
                        f_obj.write("parsed 2nd_AnnotSV output" + "\n")
                else:
                    raise subprocess.CalledProcessError()



except:
    print(
        "the status file is not found or parse AnnotSV output cannot be performed,\n ")
    exit()


