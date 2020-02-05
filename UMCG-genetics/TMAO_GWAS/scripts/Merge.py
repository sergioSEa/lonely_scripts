import sys
from pathlib import Path

Dir = sys.argv[1]

#command = "java -jar /apps/software/picard/2.9.2-Java-1.8.0_74/picard.jar GatherVcfs "
command = "bcftools concat -o {Filtered} ".format(Filtered=sys.argv[2])
for vcf in Path(Dir).glob("*.vcf"):
	if "Merged" in str(vcf): continue
	command += str(vcf) + " "
	#command += "I=" +str(vcf) + " "
#command += "O=PCA_input.vcf"



from subprocess import call
print(command)
call(command, shell=True)

