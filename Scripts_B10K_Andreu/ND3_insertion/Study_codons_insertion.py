


From = 163 
To = 180

C0 = []
C1 = []

Input = "MSA_ins.fa"
with open(Input, "r") as INPUT:
	for line in INPUT:
		if line[0] == ">" :
			codons_0 = []
			codons_1 = [] 
			continue
	
		Region = line[From-1:To]
		Region_frameshift = line[From-1+9:From+11] + line[From+12:To+1]
		
		#Check Codons in 0 reading frame
		n = 0 
		while n+3 <= len(Region):
			Codon = Region[n:n+3]
			codons_0.append(Codon)
			if n + 3 <= len(Region_frameshift):
				Codon = Region_frameshift[n:n+3]
				codons_1.append(Codon)
			n +=3
		C0.append(codons_0)
		C1.append(codons_1)

T = {}
T2 = {}
for n in range(0,len(codons_0)):
	Position = [i[n] for i in C0]
	dic_seen = {}
	Total = 0
	for P in Position:
		if P not in dic_seen: dic_seen[P] = 0
		dic_seen[P] += 1
		Total += 1
	for key in dic_seen:
		dic_seen[key] = dic_seen[key]/float(Total)

	T[n] = dic_seen

	if n < len(codons_1):
		Total = 0
		Position = [i[n] for i in C1]
		dic_seen = {}
		for P in Position:
			if P not in dic_seen: dic_seen[P] = 0
			dic_seen[P] += 1
			Total += 1
		for key in dic_seen:
			dic_seen[key] = dic_seen[key]/float(Total)

		T2[n] = dic_seen





print(T)
print(T2)
