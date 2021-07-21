File=$1
R2=$2

OUTPUT=$3

TEMP="temporal_recoding$4"
TEMP2="temporal_recoding2$4"
TEMP3="temporal_recoding3$4"



awk '/^[^#]/ {OFS = "\t"; if($14=="NA") {a==NA;} else if ($14=="0") {a=0;} else if ($14=="SQRT_V_STAT") {a="SQRT_V_STAT"} else {a=1/$14} ; print $1":"$2, $1, $2, $4, $3, "NA", $15, a, $16, $6, $7, $9, $8, $5}' <(gzip -dc $File) > $TEMP
tail -n +2 $TEMP > $TEMP2


paste --delimiters='\t' $TEMP2 $R2 > $TEMP3

cat HEADER $TEMP3 > $OUTPUT


