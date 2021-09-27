#1. Get a conda envirounment with SRA-toolkit: https://anaconda.org/bioconda/sra-tools   ; https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
/groups/umcg-tifn/tmp01/users/umcg-dwang/opt/Anaconda3/condabin/conda activate /groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/conda_env/ATLAS
#2. by default SRA toolkit saves SRA files in /home ;  as it doesnt fit we need to change the default directory: https://www.biostars.org/p/159950/
echo '/repository/user/main/public/root = "/^Cp"' > $HOME/.ncbi/user-settings.mkfg
#3. Download SRA (prefetch) into tmp, and convert to fastq
while read SRA ; do
        prefetch $SRA
        fastq-dump --split-files --gzip  /tmp/sra/$SRA\.sra
        rm /tmp/sra/$SRA\.sra
done < BGC300_ascesions.txt #Get this file from the SRA bioproject. Press "SRA Experiments" to see all SRA IDs for each sample. Download them all from "Send to" > "File" > Format:"Accession List"     

