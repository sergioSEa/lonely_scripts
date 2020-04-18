From = 163 
To = 180

C0 = []
C1 = []


Input = "MSA_ins.fa" #Output from Divide.py
with open(Input, "r") as INPUT:
	for line in INPUT:
		if line[0] == ">" :
			Bisho = line
			continue
		Region = line[From-1:To]
		Region_frameshift = line[From-1+9:From+11] + line[From+12:To+1]
		if line[175:178] == "gtg": exit(Bisho)
		continue
		#Check Codons in 0 reading frame
		n = 0 
		ITERATION = 1
		ITERATION2 = 1
		while n+3 <= len(Region):
			if len(C0) < ITERATION: C0.append({})
			Codon = Region[n:n+3]
			if "n" in Codon: 
				n += 3
				continue
			if Codon not in C0[ITERATION-1]:
				C0[ITERATION-1][Codon] = 0
			#if ITERATION-1 == 3:
			#	if Codon == "atc":
			#		print(Codon, line)
			C0[ITERATION-1][Codon] += 1
			ITERATION += 1
			if n + 3 <= len(Region_frameshift):
				if len(C1) < ITERATION2: C1.append({})
				Codon = Region_frameshift[n:n+3]
				if "n" in Codon: 
					n += 3
					continue
				if Codon == "tcc":
					exit([Codon,line])
				if Codon not in C1[ITERATION2-1]:
					C1[ITERATION2-1][Codon] = 0
				C1[ITERATION2-1][Codon] += 1
				ITERATION2 += 1
			n += 3


l = ["Sequence","Codon_position","Codon","Counts"]
print("\t".join(l))
for n in range(len(C0)):
	for item in C0[n]:
		l = ["regular",str(n+1),str(item),str(C0[n][item])]
		print("\t".join(l))
	if n < len(C1):
		for item in C1[n]:
			l = ["frameshift",str(n+1),str(item),str(C1[n][item])]
			print("\t".join(l))

