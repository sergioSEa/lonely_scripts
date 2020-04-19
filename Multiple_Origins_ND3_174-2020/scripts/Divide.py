
#Input: MSA of diapsida and table with status in position 174 in diapsida
MSA = "Alignment_diapsida.fa"
TSV = "Diapsida_table.tsv"


gaps = []
insertion = []
with open(TSV) as File:
	for line in File:
		l = line.rstrip().split()
		if l[1]  == "-":
			gaps.append(l[0])
		else:
				
			insertion.append(l[0])

GAP = ""
INS = "" 
with open(MSA) as File:
	for line in File:
		if line[0] == ">":
			l = line.rstrip()[1:]
			if l in gaps:
				add_gap = True
				add_insertion = False
			elif l in insertion:
				add_insertion = True
				add_gap = False
			else:
				print(l)
				exit("Error")
		if add_gap == True:
			GAP += line			
		elif add_insertion == True:
			INS += line



with open("MSA_gap.fa","w") as M:
	M.write(GAP)
with open("MSA_ins.fa","w") as M:
	M.write(INS)
