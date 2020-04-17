#Input: 
#taxonomy table for each record
#B10K records
#Table with with the nucleotides in position 174 for all sequences
taxonomy = "taxonomy"
B10K_data = "B10K.fa"
Final_set = "Position_nucleotide.tsv"

OUT = "Position_nucleotide2.tsv"
OUTPUT = open(OUT,"w") 
OUTPUT.close()


#B10K ID to LAtin name
translation = {}
with open("NAME-ID.txt") as File:
        for line in File:
                if "ID" in line: continue
                l = line.split(",")
                translation[l[0]] = l[1].rstrip().replace(" ","_")

#Get which records belong to b10k
list_b10K = []
with open(B10K_data) as B10K:
	for line in B10K:
		if line[0] == ">":
			l = line[1:-1].split("_")[0]
			l = translation[l]
			list_b10K.append(l)

dic_taxa = {}
# Hash of Taxonomical_name (key) and taxonomic_level (value)
with open(taxonomy) as F:
	for line in F:
		l = line.rstrip().split()
		dic_taxa[l[0]] = l[1]

# If record is the B10K ones tehn it is a bird. Check the others.
with open(Final_set) as F:
	for line in F:
		l = line.rstrip().split()
		if l[0] in list_b10K:
			taxon = "Aves"
		elif l[0] in dic_taxa:
			taxon = dic_taxa[l[0]]
		else:
			print(line)
			continue

		new_line =  "\t".join([l[0],l[1],taxon,"\n"])
		with open(OUT, "a") as OUTPUT:
			OUTPUT.write(new_line)	 

