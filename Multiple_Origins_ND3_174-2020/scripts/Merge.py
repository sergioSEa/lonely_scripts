from Bio import SeqIO

#Given two sets of ND3 (NCBI and B10K), merge them and remove duplicates
# In-house B!0K records need to change IDs by actual taxonomical name




#####Get the translation to latin name from B10K ID
translation = {}
with open("../../NAME-ID.txt") as File:
	for line in File:
		if "ID" in line: continue
		l = line.split(",")
		translation[l[0]] = l[1].rstrip().replace(" ","_")

####Save in a dic the B10K sequences. Latin_name : Sequence
Available = {}
with open("Aligned_ND3_seqs_B10K.fa") as B10K:
	S = SeqIO.parse(B10K, "fasta")
	for Species in S:
		N = Species.name.split("_")[0]
		#N0 = N
		N = translation[N]
		Available[N] = str(Species.seq).upper()




####Go through ND3 NCBI records. If they are already in B10K, add them to a list of already done. If they  have been already seen between the NCBI ones, add them to a second list of repeated in NCBI. 
#Otherwise: Attach to a string including all seqs

FINAL = ""
Repeated = []
Done_NCBI = []
repeated_NCBI = []
with open("ND3_genes.fa","r") as File:
	Entry = SeqIO.parse(File,"fasta")
	for Species in Entry:
		N = Species.name.split("-")[0]
		if N in Available:
			Repeated.append(Species.name)
			continue
		else:
			if N in Done_NCBI:
				repeated_NCBI.append(Species.name)
				continue
		Done_NCBI.append(N)
		ADD = ">{NAME}\n{sequence}\n".format(NAME=N, sequence=str(Species.seq))	
		FINAL += ADD
##Attach the B10K in the FINAL (the NCBI which passed all restrictions).
for key in Available:
	ADD = ">{NAME}\n{sequence}\n".format(NAME=key, sequence=str(Available[key]))
	FINAL += ADD


##Get some numbers

print("Number of unique records: " + str(len(Available.keys())))
print("Number of records from NCBI already in B10K:" + str(len(Repeated)))
print("Number of records from NCBI who were removed because repeated:" + str(len(repeated_NCBI)))


with open("ND3_merged.fa","w") as OUT:
	OUT.write(FINAL)
with open("Overlapping_B10K-NCBI.txt", "w") as OUT:
	OUT.write("\n".join(Repeated))
with open("Repeated_NCBI.txt", "w") as OUT:
	OUT.write("\n".join(repeated_NCBI))
