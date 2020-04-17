import sys
FILE = sys.argv[1]

save = []
with open(FILE) as MSA:
	for line in MSA:
		if line[0] == ">":
			save.append(line)
			continue
		
		if line[175-1] == "-" and line[174-1] != "-":
			line = line[0:(174-1)] +"-"+ line[174-1] + line[176-1:]
		save.append(line)
print("".join(save)) 	
