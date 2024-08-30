from pathlib import Path
import subprocess 
import argparse


def Find_reads(Dir, Suffix, FR="_R1", RR="_R2"):
    '''
    Designed to search for files within a directory (and its subdirectories) that match a specific suffix 
    '''
    if Suffix[0] != "." : Suffix = "." + Suffix
    Dic = {}
    for File in Path(Dir).rglob("*" + Suffix):
        FileName = File.name.replace(Suffix, "")
        Name = FileName.replace(FR, "").replace(RR, "")
        if Name not in Dic: Dic[Name] = {}
        if FR in FileName: Dic[Name]['Forward'] = str(File)
        elif RR in FileName: Dic[Name]['Reverse'] = str(File)
    return(Dic)

def count_reads_in_fastq(fastq_file):
    if '.gz' in str(fastq_file):
        #command = f"zgrep -c '^@' {fastq_file}"
        command =  f'zcat {fastq_file} | wc -l'
    else:
        command =  f'cat {fastq_file} | wc -l'
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # Check if there was an error
    if result.returncode != 0:
        print("An error occurred:", result.stderr)
        return None
    # Return the count as an integer
    return int(result.stdout.strip())/4
def count_reads_in_fastq_bio(fastq_gz_file):
    # Open the gzipped FASTQ file
    with gzip.open(fastq_gz_file, "rt") as handle:
        # Use SeqIO to parse the file and count the number of records
        read_count = sum(1 for _ in SeqIO.parse(handle, "fastq"))
    return read_count
        
def Get_readNumber(Dic, Method = "Biopython", Reverse=False ):
    print("Available methods: [Biopython, Count]")
    Methods = ["Biopython", "Count"]
    if Method not in Methods: return()
        
    if Method == "Biopython":
        from Bio import SeqIO
        import gzip
        Function_do = count_reads_in_fastq_bio
    elif Method == "Count":
        Function_do = count_reads_in_fastq
   
    
    Dic_number = {}
    for Name in Dic:
        Dic_number[Name] = {}
        
        F = Dic[Name]['Forward']
        R = Dic[Name]['Reverse']
        F_reads = count_reads_in_fastq(F)        
        Dic_number[Name]["Forward"] = F_reads
        print(Name, str(F_reads) )
        if Reverse == True:
            R_reads = count_reads_in_fastq(R)
            Dic_number[Name]["Reverse"] = R_reads
    
    return(Dic_number)
def Save_counts(Dic_counts, OutFile="FQ_counts.txt" ):
    with open(OutFile, "w") as O:
        O.write("Sample_ID\tCounts\tMillion_counts\n")
        for Sample in Dic_counts:
            Count = Dic_counts[Sample]["Forward"]
            Save = Sample + "\t" + str(Count) + "\t" + str(float(Count)/1000000) + "\n"
            O.write(Save)
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute read depth from FQ files.')

    parser.add_argument('--dir', '-d', required=True, help='Directory where to check for FQ files')
    parser.add_argument('--suffix', '-s', default='fastq.gz', help='Suffix for FQ files, e.g.: fq.gz or fastq.gz')
    parser.add_argument('--biopython', action='store_true', help='Compute the read number with Biopython (flag), otherwise a system call will be used')
    parser.add_argument('--forward_read', default='_R1', help='How are forward reads specified (default: "_R1")')
    parser.add_argument('--reverse_read', default='_R2', help='How are reverse reads specified (default: "_R2")')
    parser.add_argument('--out_file', default='FQ_counts.txt', help='Name of file where to save the output')
    
    args = parser.parse_args()
    FQ = Find_reads(args.dir, args.suffix, FR=args.forward_read, RR=args.reverse_read)
    if args.biopython:
        Counts = Get_readNumber(FQ, Method = "Biopython")
    else:
        Counts = Get_readNumber(FQ, Method="Count")
    Save_counts(Counts, OutFile=args.out_file)

