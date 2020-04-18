###Computes Shanon entropy per position from each lineage

import numpy as np

MSA = "Aligned_clean.fa"
MSA = "../Data/Corrected_Alignment_reptiles.fa"
METADATA = "../Data/Metadata_sp.csv"


Dic_frequencies = {"Testudines":{}, "Birds":{}, "Crocodile":{}, "Lizards":{}}

def Shanon_entropy(column_vector):
	Abecedary = ["A","T", "G", "C"]
	T = 0
	All_F = []
	for N in Abecedary:
		F = column_vector.count(N)
		T += F
		All_F.append(float(F))
	All_F = np.asarray(All_F)
	if T == 0:
		return(np.nan)
	All_F = All_F/T
	All_F[All_F==0] = 1
	logs = np.log2(All_F)
	E = np.multiply(All_F, logs)
	Entropy = -1 * sum(E)
	return(Entropy)


meta_information = {}
with open(METADATA) as M:
	for line in M:
		if "X1" in line: continue
		l = line.rstrip().split(",")
		l = [ x.strip('"') for x in l ]
		meta_information[l[1]] = l[2:]
 


append_next = False
with open(MSA) as M:
	for line in M:
		line = line.rstrip()
		if line[0] == ">":
			if line.lstrip(">") not in meta_information: continue
			M = meta_information[line.lstrip(">")]
			
			animal = M[-2]
			Gap = M[-1]
						
			if Gap not in Dic_frequencies[animal]: Dic_frequencies[animal][Gap] = []
			append_next = True
		else:
			if append_next != True: continue
			Dic_frequencies[animal][Gap].append(line)
			append_next = False


Outcome = ["\t".join(["Position", "Entropy", "Group", "Gap_status", "codon_bin"])]

for Animal in Dic_frequencies:
	for Gap_status in Dic_frequencies[Animal]:
		prev_entropy = []
		for n in range(len( Dic_frequencies[Animal][Gap_status][0])):
			Position =  [ x[n].upper() for x in Dic_frequencies[Animal][Gap_status] ]
			Entropy = Shanon_entropy(Position)
			if "-" in str(Entropy): Entropy = 0
			
			if (n+1)%3==0:
				prev_entropy.append(Entropy)
				if np.nan not in prev_entropy: codon_entropy = sum(prev_entropy) / 3
				else:
					codon_entropy = 0
					m = 0
					for i in prev_entropy:
						if str(i) != "nan":
							codon_entropy += i
							m += 1
					try:
						codon_entropy = codon_entropy / m
					except:
						codon_entropy = np.nan
					
				prev_entropy = []
			else:
				prev_entropy.append(Entropy)
				codon_entropy = np.nan
			Outcome.append("\t".join([str(n+1), str(Entropy), Animal, Gap_status, str(codon_entropy)]))

OUTPUT = "\n".join(Outcome)

print(OUTPUT)