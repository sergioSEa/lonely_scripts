
File_with_header=$1
input_file=$2

mkdir temp
head -n 1 $File_with_header > temp/temp_header

touch temp/temp_EMPTY
paste -d "\t" temp/temp_EMPTY temp/temp_header > temp/header


awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' temp/header > temp/transposed_header

printf "\n" > temp/temp_EMPTY
cat temp/temp_EMPTY temp/transposed_header > temp/temp_transposed

cat temp/header $input_file > temp/temp_output


paste -d "\t" temp/temp_transposed temp/temp_output > $input_file

rm -r temp

