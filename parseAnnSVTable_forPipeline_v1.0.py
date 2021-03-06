#######################################################
# Author :  Jainy Thomas
# date   :  June 2019
# email  :  jainythomas1@gmail.com
# Purpose :  parse AnnotSV output for sorting and printing variants inherited in interested individuals, and variants of particular interest
#
#####################################################
changelog = '''changelog:
  - v1.0 = 1 June 2019

            This script is used for formating and modifying commandline options
  - v2.0  = 4 June 2019
            get the column and details from commandline
            hardcoded for type of variant picked (e.g not in intron)
            also hardcoded for individuals of 12, if there is change from the script needs to be modified
            also SVinfo is not included in the ANNOTSV table, if that is included the script needs to be modified to 
            account for the changes in the column numbers
  - v3.0 = June 25 2019
            previous version was for especially for choosing introns, This version I am trying to choose that annotation to select as an optional step
  - v4.0 = February 11 2020
            modify to parseANNOTSV outputput from the TEanalysis pipeline
            
            # Identify the candidates and print whole lines to a file
            # also print file which can be used as input file for validation script in perl  
            
                  

'''
import argparse
import re

###########################_LOAD AND CHECK_#########################################
####################################################################################
parser = argparse.ArgumentParser()
#parser.add_argument("-v","--verbosity",help="increase verbosity",action="store_true")
parser.add_argument("-t","--table",help="file containing the ANNOTSV output onlyfull")
#parser.add_argument("-t","--table",type=argparse.FileType('r'),help="file containing the ANNOTSV output")
parser.add_argument("-i",'--individuals',help=" column of individuals that are need to be share the genotype either "
                                              "as het or homo starting first column of the individuals as 1;"
                                              "(./. is considered as absent)",type=str)  #e.g, 3.4.5 (here in the sample for individuals 122751,1121204,120201)
parser.add_argument("-a",'--noshare',help="not sharing (./. also considered as none",type=str,nargs='?',default=0) # 9 for the family containing 122583
parser.add_argument("-n","--total",help="total number of individuals genotyped",type=int,nargs='?',default=12)
parser.add_argument("-j","--startcol",help="starting column number for individuals",type=int,nargs='?',default=9)
# add the option on parse based on scores
parser.add_argument("-s","--scores", help = " parse on scores, select those with scores 3 or higher", type=int, nargs='?', default=3)
parser.add_argument("-o","--outputfile", help = "outputput file")
parser.add_argument("-vhg", "--hg_ver", help="version of the human genome")

#get column numbers of individuals or default coli==9
#get total number of individuals nind=12
args= parser.parse_args()

if not args.outputfile:
    args.outputfile = args.table
#numlist = args.individuals.split(',')
if args.individuals != "None":
    sharelist = [int(item) for item in args.individuals.split(',')]
    nofinds= len(sharelist)
    print(sharelist, "shared individuals")
print("The args noshare is:", args.noshare)
print("The args total is:",args.total)
if (args.noshare != 0):
    no_sharelist = [int(item) for item in args.noshare.split(',')]
    print(no_sharelist, "non-shared individuals")
    noabsind=len(no_sharelist)



#print(args.table)
######################################_FUNCTIONS_##########################################
###########################################################################################

def printofile_dictnary(dictn,filename,header,outtype):
    outfilobj = open(filename + "_" + outtype + ".out", "w")
    outfilobj.write(header + '\n')
    for elements in dictn:
        svline = dictn[elements]
        #print(svline)
        outfilobj.write(svline + '\n')
    outfilobj.close()

def check_gt_similiarity(list,dict):
    #totnumid=len(list)
    z=0
    for k in list:
        if dict[k] == "0/1" or dict[k] == "1/1":
            z += 1
        else:
            break
    return(z)


def check_gt_difference(list,dict):
    totnumid=len(list)
    #print(list)
    z=0
    for k in list:
        if dict[k] == "0/0" or dict[k] == "./.":
            z += 1
        else:
            break
    return(z)

def print_keydict(dictn,filename,header,outtype):
    outfilobj = open(filename + "_" + outtype + ".out", "w")
    outfilobj.write(header + '\n')
    for elements in dictn:
        #svline = dictn[elements]
        #print(svline)
        outfilobj.write(elements + '\n')
    outfilobj.close()



##########################################_MAIN_###########################################
###########################################################################################


herit={}
highscorecand = {}
hscore_candvali = {}
tsdcheckcandi = {}
#linenumber = 0
with open(args.table) as file_obj:
    global nucols
    startcolu = (args.startcol - 1)
    lastcolu = startcolu + args.total
    for line in file_obj:#read line by line
        #linenumber +=1
        line = line.rstrip()#remove newline character
        found = re.search(r"AnnotSV ID", line)
        cols = line.split("\t")

        if found:
            #print(found)
            title = line
            dgt = {}
            nucols = len(cols)
            #print(nucols)
            h=1
            for k in range(startcolu,
                           lastcolu):  # changed to automatically detect the number of individuals that are genotyped
                # print(i)
                keyind = "indi_" + str(h)

                valueindgt = cols[k]
                dgt[keyind] = valueindgt

                h = h + 1
                k = k + 1


        else:
            chr = cols[1]
            start = cols[2]
            end = cols[3]
            uniqueid = cols[0]
            if args.individuals != "None":

                gt = {}

                # icolnumber = 9

                n = 1

                for i in range(startcolu,
                               lastcolu):  # changed to automatically detect the number of individuals that are genotyped
                    # print(i)
                    #keyind = "indi_" + str(n) + "_gt"

                    valueindgt = cols[i].split(":", 1)[0]
                    #dgt[keyind] = valueindgt
                    gt[n] = valueindgt
                    n = n + 1
                    i = i + 1
                #print(dgt)
                #print(gt)
                if (args.noshare != 0):
                    #print("Here I am", args.noshare)
                    nofmatches = check_gt_similiarity(sharelist, gt)
                    #print(nofmatches)
                    nofabsgt = check_gt_difference(no_sharelist, gt)# checking to see if any one of the individuals
                    # identified in the command line is not to share is having the presence
                    if (int(nofmatches) == int(nofinds)) and int(nofabsgt) == int(noabsind):
                        herit[uniqueid] = line
                else:
                    nofmatches = check_gt_similiarity(sharelist, gt)
                    if (int(nofmatches) == int(nofinds)):
                        # if ((indi_3_gt == "0/1" or indi_3_gt == "1/1") and (
                        # indi_4_gt == "0/1" or indi_4_gt == "1/1") and (
                        # indi_5_gt == "0/1" or indi_5_gt == "1/1")):
                        herit[uniqueid] = line
                    else:
                        continue
            if (args.scores):
                scorecol = int(nucols)-1
                scoreval = (cols[scorecol])#    continue
                if (int(scoreval) >= args.scores):
                    highscorecand[uniqueid] = line

                    if args.hg_ver == "hg38":
                        modchr = "chr" + chr
                        tsdkey = "\t".join([modchr, start, end, cols[5]])
                        tsdcheckcandi[tsdkey] = 1
                    elif args.hg_ver == "hg19":
                        tsdkey = "\t".join([chr, start, end, cols[5]])
                        tsdcheckcandi[tsdkey] = 1
                    for inds in dgt:
                        inname = dgt[inds]
                        newid = uniqueid + "_" +inname
                        newvalue = "\t".join([inname,chr,start,end,cols[5]])
                        hscore_candvali[newid]=newvalue


                        #going to make a table to be used by validation script currently designed to do for all individuals
                    # for loop for candidate repeat for all individuals



file_obj.close()
if (args.individuals):
    printofile_dictnary(herit,args.outputfile,title,"segre")
if (args.scores):
    printofile_dictnary(highscorecand,args.outputfile,title,"score")
    inputtabletitle = "#indivbam\tchr\tstart\tend\tTE"
    printofile_dictnary(hscore_candvali, args.outputfile, inputtabletitle, "input_table")
    title2 = "#chr\tstart\tend\tTE"
    print_keydict(tsdcheckcandi, args.outputfile, title2, "for_TSDeval")

