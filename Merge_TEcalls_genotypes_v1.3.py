#######################################################
# Author :  Jainy Thomas
# date   :  June 2019
# email  :  jainythomas1@gmail.com
# Purpose :  Merge SV list with genotype in bed format
#
#####################################################
changelog = '''changelog:
  - v1.0 = 5 June 2019
			Input file format 
			1	1594397	1660860	DUP	N	<DUP>	GT:SU:PE:SR:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB	0/1:10:10:0:14:14.52:-5,-4,-11:121:111:9:111:9:0:0:0:111:9:0.075	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	0/1:15:15:0:45:92.38:-11,-2,-6:112:97:14:97:14:0:0:0:97:14:0.13	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.	./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.
  - v 2.2 = 6 December 2019
  			Added to option get the total number of columns from the bed file
  			The script can now accept the file with headers with a pound sign (#) at the beginnning and print the headers in the output
  - v 2.3 = 9 December 2019
  			can use with file containing duplicate lines
  			fixed a bug when printing last line of the file when line ends in line number 3
  -v 2.4	10 December 2019
  			sort the dictionary to print in the correct order	using pandas
			To do:
			#I could do store the data in multidimensional array by using numpy, will get to that later
  -v 3.0    15 January 2020
  			Modify the script to adapt the output from MELT and Transurveyor
  -v 1.2    June 3 2020
  			remove the criteria that type of TE identfied should be same.
  			intrduced the ability to alter the length of the bp to merged
  			takes the MELT genotype if available
  -v 1.3    Sep 10 2020	
  			Not to merge, if the two different of insertions are nearby especially identified by MELT	
			
'''
usage ='''


	   '''

import argparse
import re
import subprocess
import pandas as pd


####################################_LOAD AND CHECK_#######################################
###########################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("-t","--table",help="file containing SVs to merge (sorted bed format)")
parser.add_argument("-l","--len",help="length of bp to merged",type=int,nargs='?',default=100)
parser.add_argument("-c","--numcol", help = "total number of columns in the bed file",type =int)
parser.add_argument("-o","--outputfile", help = "outputput file")
args= parser.parse_args()

#global args.numcol

#print(args.table)
if not args.outputfile:
	args.outputfile = args.table + ".merged.sorted.bed"
if not args.len:
	args.len = 100
######################################_FUNCTIONS_##########################################
###########################################################################################

def comparetwolines(llist1, llist2, lline1, lline2, lintyp,lenbp):

	#print(llist1, llist2, lintyp,lenbp)
	#print("secondstation")



										#start2   -          #start1                         #end2    -    #end1                          #end1    -            end2
	if ((llist1[0] == llist2[0]) and (abs(int(llist2[1]) - int(llist1[1])) <= lenbp) and ((abs((int(llist2[2]) - int(llist1[2]))) <= lenbp) or (abs((int(llist1[2]) - int(llist2[2]))) <=lenbp))):
		# print(" #start2   -    #start1 ")
		# print(abs(int(llist2[1]) - (int(llist1[1])) <= 100))
		# print(" #end2    -    #end1 ")
		# print(((int(llist2[2]) - int(llist1[2])) <= 100))
		#
		# print(" #end1   -  #end2 ")
		# print(((int(llist1[2]) - int(llist2[2])) <= 100))
		# print(abs(int(llist1[2])-int(llist2[2])))
		cols1 = lline1.split("\t")
		cols2 = lline2.split("\t")
		if cols1[5].isalpha() and "M" in cols1[4] and cols2[5].isalpha() and "M" in cols2[4] and cols1[5] != cols2[5] and cols1[6] != cols2[6]:
			# print("fourth station")
			if lintyp == "formatted":
				return ("unique", lline2)
			elif lintyp == "new":
				return ("unique", lline1)

		else:
			startf = llist1[1]
			if llist2[2] > llist1[2]:
				endf = llist2[2]
			elif llist2[2] == llist1[2]:
				endf = llist1[2]
			else:
				endf = llist1[2]

			# type of tool: here M or T
			tool1 = cols1[4]
			tool2 = cols2[4]
			if tool1 == tool2:
				toolf = tool1
			elif tool1 in tool2:
				toolf = tool2
			elif tool2 in tool1:
				toolf = tool1
			else:
				toolf = tool1 + ";" + tool2

			# 5 is TSD in melt score in TR
			iinfo1 = cols1[5]
			iinfo2 = cols2[5]
			print(iinfo1, iinfo2)
			# cont_digit1 = any(map(str.isdigit,iinfo1))
			# cont_digit2 = any(map(str.isdigit, iinfo2))
			if iinfo1 == iinfo2:
				iinfof = iinfo1
			elif iinfo1.isalpha() and iinfo2.isalpha():
				if iinfo1 == "null":
					iinfof = iinfo2
				elif iinfo2 == "null":
					iinfof = iinfo1
				elif iinfo1 != "null" and iinfo2 != "null":
					iinfof = iinfo1 + "second" + iinfo2
			elif iinfo1.isalpha() and any(map(str.isdigit, iinfo2)) and iinfo1 != "null":
				iinfof = iinfo1
			elif iinfo1.isalpha() and any(map(str.isdigit, iinfo2)) and iinfo1 == "null":
				iinfof = iinfo2


			elif any(map(str.isdigit, iinfo1)) and iinfo2.isalpha() and iinfo2 != "null":
				iinfof = iinfo2
			elif any(map(str.isdigit, iinfo1)) and iinfo2.isalpha() and iinfo2 == "null":
				iinfof = iinfo1
			elif any(map(str.isdigit, iinfo1)) and any(map(str.isdigit, iinfo2)):
				iinfof = iinfo1  # picking the firstone #

			#
			# if re.search("_",cols1[5]):
			# 	discore1 = cols1[5].split("_")
			# else:
			# 	tsd = cols1[5]
			# if re.search("_",cols2[5]):
			# 	discore2 = cols2[5].split("_")
			# else:
			# 	tsd = cols2[5]
			# 6 is strand in melt
			sinfo1 = cols1[6]
			sinfo2 = cols2[6]
			if sinfo1 == "+" and sinfo2 == "+":
				sinfof = sinfo1
			elif sinfo1 == "-" and sinfo2 == "-":
				sinfof = sinfo2
			elif sinfo1 == "+" or sinfo1 == "-" and sinfo2 != "+" or sinfo2 != "-":
				sinfof = sinfo1
			elif sinfo2 == "+" or sinfo2 == "-" and sinfo1 != "+" or sinfo1 != "-":
				sinfof = sinfo2
			elif sinfo1 != "+" or sinfo1 != "-" and sinfo2 != "+" or sinfo2 != "-":
				sinfof = sinfo1

			# 7 is Meinfo
			Meinfo1 = cols1[7]
			Meinfo2 = cols2[7]
			if Meinfo1 == "NA" and Meinfo2 == "NA":
				Meinfof = "NA"
			elif Meinfo1 == "NA" and Meinfo2 != "NA":
				Meinfof = Meinfo2
			elif Meinfo2 == "NA" and Meinfo1 != "NA":
				Meinfof = Meinfo1
			elif Meinfo1 != "NA" and Meinfo2 != "NA":
				if Meinfo1 in Meinfo2:
					Meinfof = Meinfo2
				elif Meinfo2 in Meinfo1:
					Meinfof = Meinfo1

				else:
					Meinfof = Meinfo1 + ";" + Meinfo2

			colsf = [llist1[0], startf, endf, llist1[3], toolf, iinfof, sinfof,
					 Meinfof]
			# print(len(cols1))
			# print(len(cols2))
			#     uniqueidf = chr1 + startf
			# 	print("third station")
			# print ("The lines comparing are\n", lline1, "\n", lline2)
			# for i in range(7, 18):
			i = 8
			while i < int(args.numcol):
				# print(i)

				nc1 = re.search(r"^\.\/\..*", cols1[i])
				nc2 = re.search(r"^\.\/\..*", cols2[i])
				gt1 = re.search(r"(^0\/.*|^1\/.*)", cols1[i])
				gt2 = re.search(r"(^0\/.*|^1\/.*)", cols2[i])
				# print("gt2 is", gt2.group(1), "nc1 is", nc1.group(1), "gt1 is", gt1.group(1), "nc2 is", nc2.group(1))

				if (gt1) and (gt2) and (gt1.group(1) == gt2.group(1)):
					# print ("The two genotypes are equal: gt1 is",gt1,"gt2 is",gt2)
					colsf.append(cols1[i])
				elif (gt1) and (gt2) and gt1.group(1) != gt2.group(1):

					# print("The two genotypes are not equal: gt1 is", gt1, "gt2 is", gt2)
					gtdata1 = gt1.group(1)
					gtdata2 = gt2.group(1)
					#print(gtdata1, gtdata2)
					predgt1 = gt1.group(1).split(":")
					predgt2 = gt2.group(1).split(":")
					if len(predgt1) == 2 and len(predgt2) == 1:
						colsf.append(cols1[i])
					elif len(predgt2) == 2 and len(predgt1) == 1:
						colsf.append(cols2[i])
					elif len(predgt2) == len(predgt1):
						#print(predgt1[0], predgt2[0])
						if predgt1[0] == predgt2[0]:

							# cols1[i] = cols1[i] + "_twogt"

							colsf.append(cols1[i])
						elif (predgt1[0] == "0/0" and predgt2[0] == "0/1") or (
								predgt1[0] == "0/1" and predgt2[0] == "1/1") or (
								predgt1[0] == "0/0" and predgt2[0] == "1/1"):
							cols2[i] = cols2[i] + "_twogt"
							colsf.append(cols2[i])
						elif (predgt1[0] == "0/1" and predgt2[0] == "0/0") or (
								predgt1[0] == "1/1" and predgt2[0] == "0/1") or (
								predgt1[0] == "1/1" and predgt2[0] == "0/0"):
							cols1[i] = cols1[i] + "_twogt"
							colsf.append(cols1[i])

						global nuofconflicts
						nuofconflicts += 1
				elif ((gt1) and (nc2)):
					colsf.append(cols1[i])
				elif ((gt2) and (nc1)):
					# print("gt2 is", gt2, "nc1 is", nc1)
					colsf.append(cols2[i])
				elif (nc2) and (nc1):
					colsf.append(cols2[i])
				#
				i += 1
			linef = "\t".join(colsf)  #
			#print(linef)
			#
			#     # compare f to line3
			return ("merged", linef)



	else:
		#print("fourth station")
		if lintyp == "formatted":
			return ("unique", lline2)
		elif lintyp == "new":
			return ("unique",lline1)

def printofile_dictnary(dictn,filename,header):
	outfilobj = open(filename + ".out.bed", "w")
	outfilobj.write(header + '\n')
	for keyelements in sorted(dictn.keys()):
		svline = dictn[keyelements]
		#print(svline)
		outfilobj.write(svline + '\n')
	outfilobj.close()

def convert_dict_to_dataframe_print(dict,filename,header):
	#outfilobj = open(filename + ".out.sorted.bed", "w")
	outfilobj = open(filename, "w")
	#print("create data from dictionary,\n")
	dfObj = pd.DataFrame.from_dict(dict, orient='index',columns =["SVLine"])# creating using dict_Keys as index

	newcol = dfObj["SVLine"].str.split("\t", expand = True,)
	#print(header)
	colnames = header.split("\t")
	#print(colnames)
	newcol.columns = colnames

	# print(newcol)
	# print(newcol.columns)
	#newcol.drop(columns = 0)
	#print(newcol.dtypes)
	newcol2 = newcol.apply(pd.to_numeric, errors='coerce').fillna(newcol)
	#print(newcol2.dtypes)
	newcol2 = newcol2.sort_values(by=[colnames[0], colnames[1]])
	# sortedindex = newcol.sort_index()
	#print(newcol2)
	newcol2[colnames[0]] = newcol2[colnames[0]].astype(str)
	newcol2[colnames[0]] = newcol2[colnames[0]].str.replace('.0','',regex = False)

	newcol2 = newcol2.sort_values(by=[colnames[0], colnames[1]])
	newcol2.to_csv(outfilobj, sep='\t',index=False)
	outfilobj.close()
	#return filename



#########################################_MAIN_############################################
###########################################################################################
#list=[]
svnonmerged = {}
svmerged = {}
linenumber = 0
storelines ={}
lastline = subprocess.check_output('tail -n 1 {}'.format(args.table),shell = True)
lastline = lastline.decode('utf-8')
#print(lastline)
colsp = lastline.split("\t")
lastlinelist = [colsp[0], colsp[1], colsp[2], colsp[3]]
nuofconflicts = 0
#print("The last line identified is ", lastlinelist)
global title
#title ="#SV_chrom\tSV_start\tSV_end\tSV_type\tREF\tALT\tFORMAT\t120180\t120201\t121204\t122026\t122577\t122583\t122584\t122585\t122586\t122587\t122588\t122751"
with open(args.table) as file_obj:

	for line in file_obj:#read line by line
		line = line.rstrip()  #
		if line.startswith("#"):

			title = line #type str

		else:
			cols = line.split("\t")
			if linenumber == -1:
				#print("station five")
				line1 = storelines["thirdline"]
				#print("third line is new first line",line1)
				storelines["firstline"] = line1
				newcols = line1.split("\t")
				linelist1 =[newcols[0],newcols[1],newcols[2],newcols[3]]
				linenumber =2
			else:
				linenumber +=1
			if linenumber == 1:
				# print("I am line number one")
				# if line == lastline:
				# 	print("I am the last line and am line2 now, \n")
				chr1 = cols[0]
				start1 = cols[1]
				end1 = cols[2]
				svtype1 = cols[3]
				linelist1 = [chr1, start1, end1, svtype1]
				#uniqueid0=chr0+start0+end0+type0
				#uniqueid1 = chr1 + start1
				line1= line
				storelines["firstline"] = line1
			if linenumber == 2:
				# if line == lastline:
				# 	print("I am the last line and am line2 now, \n")
				chr2 = cols[0]
				start2 = cols[1]
				end2 = cols[2]
				svtype2 = cols[3]
				linelist2 = [chr2, start2, end2, svtype2]
				line2= line
				storelines["secline"] = line2
				# print("I am on line2")
				# print("The secondline is ",line2)
				#print("The last line is",lastline)
				#if line2 == lastline:#
				if linelist2==lastlinelist:
					# print("station seven")
					# print("I am comparing line1 and line 2 now, but the  line2 is lastline")
					(result, modline) = comparetwolines(linelist1, linelist2, line1, line2, "new",args.len)
					if result == "merged":
						data = modline.split("\t")
						uniqueidf = data[0] + "_" + data[1]
						#linelistf = [data[0], data[1], data[2], data[3]]
						#print(modline)
						svmerged[uniqueidf] = modline
						# print("For the last line, Line 1 and Line2 are merged now for ")
					elif result == "unique":
						dataline = modline.split("\t")
						uniqueidd = dataline[0] + "_" + dataline[1]
						#here line1 is added to nonmerged list, line2 is added outside the loop
						if (uniqueidd in svnonmerged.keys()):

							existline = svnonmerged[uniqueidd]
							currentdata = existline.split("\t")
							if (abs(int(dataline[2]) - int(currentdata[2])) > 100):
								uniqueidd = dataline[0] + "_" + dataline[1] + "_" + dataline[2]
								svnonmerged[uniqueidd] = modline
								# print("for last line Line1 is added to nonmerged now")
						else:
							svnonmerged[uniqueidd] = modline
							# print("for last line Line1 is added to nonmerged now")
				#uniqueid1 = chr1 + start1 + end1 + type1
				#uniqueid2 = chr2 + start2
			if linenumber == 3 :#if the last line, it is added outside the loop
				# if line == lastline:
				# 	print("I am the last line and am line3 now, \n")
				chr3 = cols[0]
				start3 = cols[1]
				end3 = cols[2]
				svtype3 = cols[3]
				linelist3 = [chr3, start3, end3, svtype3]
				line3 = line
				storelines["thirdline"] = line3
				# uniqueid1 = chr1 + start1 + end1 + type1
				#uniqueid3 = chr3 + start3
				# print("I am comparing line 1 and line 2 now")
				#compare two lines 1 and 2
				(result, modline) = comparetwolines(linelist1, linelist2, line1, line2, "new",args.len)
				if result == "merged":
					data = modline.split("\t")
					uniqueidf = data[0] + "_" + data[1]
					linelistf = [data[0], data[1], data[2], data[3]]
					# print(modline)
					svmerged[uniqueidf] = modline
					linef = modline
					# print("first and second lines merged now")
					#print("I am comparing merged line and line 3 now", linef,"\n",line3)
					(result, modline) = comparetwolines(linelistf, linelist3, linef, line3, "formatted",args.len)
					#print(linelist3)


					#print(result, modline)
					#print(line3)
					if result == "merged":
						data = modline.split("\t")
						uniqueidf = data[0] + "_"+ data[1]
						linelistf = [data[0], data[1], data[2], data[3]]
						svmerged[uniqueidf] = modline
						# print("second and third merged now")
						# print (modline)
						storelines["thirdline"] = modline
						linenumber = -1
					elif result == "unique":
						dataline = modline.split("\t")
						uniqueidd = dataline[0] + "_" + dataline[1]

						# print("second and third are staying unique")
						# print(modline)
						linenumber = -1

				elif result == "unique":
					dataline = modline.split("\t")
					uniqueidd = dataline[0] + "_" + dataline[1]
					if uniqueidd in svmerged.keys():
						currentmergedata = svmerged[uniqueidd]
						splitmerdat = currentmergedata.split("\t")
						if (abs(int(dataline[2]) - int(splitmerdat[2])) == 0):
							pass
						else:
							svnonmerged[uniqueidd] = modline

					elif (uniqueidd in svnonmerged.keys()):

						existline = svnonmerged[uniqueidd]
						currentdata = existline.split("\t")
						if (abs(int(dataline[2]) - int(currentdata[2])) > 100):
							uniqueidd = dataline[0] + "_" + dataline[1] + "_" + dataline[2]
							svnonmerged[uniqueidd] = modline
							#print("first and second lines are unique and line 1 is added to nonmerged")
						else:
							svnonmerged[uniqueidd] = modline
							#print("first and second lines are unique and line 1 is added to nonmerged_2")
					else:

						svnonmerged[uniqueidd] = modline
						#print("first and second lines are unique and line 1 is added to nonmerged_3")
					#print("I am comparing  line2 and line 3 now")
					#compare line 2 to line 3
					(result, modline) = comparetwolines(linelist2, linelist3, line2, line3, "new",args.len)
					if result == "merged":
						data = modline.split("\t")
						uniqueidf = data[0] + "_" + data[1]
						linelistf = [data[0], data[1], data[2], data[3]]
						svmerged[uniqueidf] = modline
						# print("second and third lines merged now but line1 was unique")
						# print(modline)
						storelines["thirdline"] = modline
						linenumber = -1
					elif result == "unique":
						#print("case unique4")
						dataline = modline.split("\t")
						uniqueidd = dataline[0] + "_" + dataline[1]
						if (uniqueidd in svnonmerged.keys()):

							existline = svnonmerged[uniqueidd]
							currentdata = existline.split("\t")
							if (abs(int(dataline[2]) - int(currentdata[2])) > 100):
								uniqueidd = dataline[0] + "_" + dataline[1] + "_" + dataline[2]
								svnonmerged[uniqueidd] = modline
								#print("second is added to nonmerged and third lines are unique and line1 was also unique")
							else:
								svnonmerged[uniqueidd] = modline
								#print("second  is added to nonmergedand third lines are unique and line1 was also unique_2")
						else:
							svnonmerged[uniqueidd] = modline
							#print("second is added to nonmerged and line1 was also unique_3")
						#print(modline)

						linenumber = -1




	# print("The last line number is ",linenumber)                #line 3 remains to compared to four
	# print("The last line in the file is ",line)
	# print("The lastoutput from the script is \n",modline)
	# print("The last result is")
	# print(result)


#it gets here only if the last line is

	if (linenumber == -1 or linenumber == 2) and result == "unique" and line == lastline:

		#print(line2)

		#print(line3)
		if linenumber == 2:
			finalline = line2
		elif linenumber == -1:
			finalline = line3
		linedata = finalline.split("\t")
		uniqueidl = linedata[0] + "_" + linedata[1]
		if (uniqueidl in svnonmerged.keys()):

			existingline = svnonmerged[uniqueidl]
			#print( "The existing line",existingline,"\n")
			datacurren = existingline.split("\t")
			if (abs(int(linedata[2]) - int(datacurren[2])) > 100):
				uniqueidl = linedata[0] + "_" + linedata[1] + "_" + linedata[2]
				svnonmerged[uniqueidl] = finalline
			else:
				svnonmerged[uniqueidl] = finalline
		else:
			svnonmerged[uniqueidl] = finalline
	#line ends unique and in line number 2, it is compared as special case for line 2-done
	#line end in unique and line number -1 (line3 already compared, so just have to be ), only one line left so can print line3 outside the loop-done
	#line end in merged and line number -1, there is nothing to be done- done
	# line ends merged and in line number 2, you might wanna compare the last line1 and line2 outside the loop may be also with the previous- done


print("The number of conflicted genotypes are",nuofconflicts)

svnonmerged.update(svmerged)#merge two dictionaries


print("Now printing the final merged file")

convert_dict_to_dataframe_print(svnonmerged,args.outputfile,title)

print('*'*50 + 'DONE'+'*'*50)

