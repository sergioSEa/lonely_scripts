from Bio import SeqIO

#Given two sets of ND3, merge them and remove duplicates
#ND3_complete_genes: Set of refseq genes
#"ND3_genes_b.fa": set of genbank genes
#Priority: Refseq > Genbank

FINAL = "" #Fasta records


repeated_NCBI = [] #If more than one record of the same species exist in the database
Available = {} 
with open("ND3_complete_genes.fa") as RefSeq:
	S = SeqIO.parse(RefSeq, "fasta")
	for Species in S:
		N = Species.name
		if N in Available:
			repeated_NCBI.append(Species.name)
		Available[N] = str(Species.seq).upper()

Repeated = []
Done_NCBI = []

##Collection of Refseq records in fasta (produced by Prepare_fasta.py
with open("ND3_genes_b.fa","r") as File:
	Entry = SeqIO.parse(File,"fasta")
	for Species in Entry:
		N = Species.name
		#Save records that are also present in the RefSeq
		if N in Available:
			Repeated.append(Species.name)
			continue
		else:
			if N in Done_NCBI: #If record not in RefSeq but already seen in Genbank
				repeated_NCBI.append(Species.name)
				continue
		Done_NCBI.append(N)
		ADD = ">{NAME}\n{sequence}\n".format(NAME=N, sequence=str(Species.seq))	
		FINAL += ADD

##Attach the RefSeq in the FINAL (the Genbank records which passed all restrictions).
for key in Available:
	ADD = ">{NAME}\n{sequence}\n".format(NAME=key, sequence=str(Available[key]))
	FINAL += ADD



with open("ND3_genes.fa","w") as OUT:
	OUT.write(FINAL)
with open("Overlapping_Refseq-Genebank.txt", "w") as OUT:
	OUT.write("\n".join(Repeated))
with open("Repeated_NCBI.txt", "w") as OUT:
	OUT.write("\n".join(repeated_NCBI))
