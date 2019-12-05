import shlex
import sys

evo = sys.argv[1]
#evo = "Passeriformes_tip_rate.tsv"
Traits = "../Analysis_Mtc_number/info/Trait_data2.csv"
translation = "../Analysis_Mtc_number/info/NAME-ID.txt"
phylo = "../Analysis_Mtc_number/info/phylo_table"


dic_trans = {}
with open(translation) as F:
	for line in F:
		if line.startswith("ID"): continue
		line = line.rstrip().split(",")
		dic_trans[line[0]] = line[1]

dic_quan = {}
n = 0
with open(evo,"r") as Q:
	for line in Q:
		line = line.rstrip()
		if n == 0:
			header = line
			n =1
			continue
		ID = line.split()[2]
		if ID == "NA": continue
		ID = ID.replace("_","-")
		Name = dic_trans[ID].replace(" ","_")
		line = "\t".join(line.split()[0:2]) + "\t" + line.split()[2].replace("_","-") + "\t" + "\t".join(line.split()[3:]) 
		
		if Name in dic_quan:
			dic_quan[Name].append(line)
		else:
			dic_quan[Name] = [line]
		

dic_p = {}
with open(phylo) as P:
	for line in P:
		line = line.rstrip()
		pl = line.split("|")
		line = "|".join(pl[0:2] + pl[48:])
		#line = "|".join(pl[0:2] + pl[20:])
		
		#ID = line.split("|")[11].replace(" ","_")
		try: ID = line.split("|")[13].replace(" ","_")
		except: exit(line)
		if line.startswith("Source"):
			header3 = line.replace("|","\t")
			continue

		line = line.replace("\t"," ")
		dic_p[ID] = line.replace("|","\t")
		


list_write = []
with open(Traits) as F:
	for line in F:
		line = line.rstrip()
		if line[-1] == ";": 
			line = line[0:-1]
		line = line.replace(";","\t")
		if line.startswith("Latin_name"): 
			Header2 = line
			continue
		#else:
		#	line = line.rstrip("\t")
		#WHITE ROWS
		try :Name = line.split()[0]
		except: continue
		
		#Birds not analyzed
		try:
			info = dic_quan[Name]
		except: continue
		
		more_info = dic_p[Name]
		
		for PROTEIN in info:
			FINAL = PROTEIN + "\t" + line + "\t"+more_info  +  "\n"
			list_write.append(FINAL)



#header = header.replace("\t","\t")
#Header2 = Header2.replace("\t","\t")
#header3 = header3.replace(";","\t")


output = sys.argv[2]
with open(output,"w") as O:
	H= header + "\t" + Header2 +"\t" + header3 +"\n"
	O.write(H)
	for line in list_write:
		O.write(line)
	#	if len(line.split("\t")) != 198:
	#		print(len(line.split("\t")))



