#!/bin/bash

###############
## HELP
#################
usage()
{
    echo "usage: make_SNPtable_from_vcf.sh [-v vcf.gz] [-c chroms] [-f firstSampleCol=10] [-m minNofCalls=2] [ -o outfileBase=<vcf> ] [ -t threads] [ -h help ]
	** note that the vcf must be gzipped **
"
}
if [ $# -lt 1 ]; then usage; exit; fi

##### Main

firstSampleCol=10
threads=1

while [ "$1" != "" ]; do
    case $1 in
        -v | --vcf )            shift
                                vcf=($1)
                                ;;
        -c | --chroms )         shift
                                chroms=($(echo ${1} | tr ',' ' '))
                                ;;
        -f | --firstSampleCol ) shift
                                firstSampleCol=$1
                                ;;
        -m | --minCalls )       shift
                                minCalls=$1
                                ;;
        -o | --outfileBase )    shift
                                outfileBase=$1
                                ;;
        -t | --threads )        shift
                                threads=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

###################
## MAIN
##################
if [ -z "$outfileBase" ]; then outfileBase=$(echo $vcf | sed 's/.vcf.gz//'); fi
echo "making snptable for chroms:[$(echo ${chroms[*]})] from $vcf, starting from column $firstSampleCol, and writing to $outfileBase.<chrom>.snptable"

makeSNPtable() {
	vcf=${1}
	firstSampleCol=${2}
	minCalls=${3}
	outfileBase=${4}
	chrom=${5}
	
	outFile=$outfileBase.$chrom.snptable; echo "working on $outFile ..."; 
	
	##PRINT HEADER
	zcat $vcf | head -1000  | grep "#C" | head -1 | awk -v fsc=$firstSampleCol '{for(i=fsc;i<=NF;i++){printf $i","}}' |  sed "s/^/${chrom},Ref,/" | sed 's/,$/\n/'  > $outFile 
	
	zcat $vcf | grep -P '^'$chrom'\t' | awk -v fsc="$firstSampleCol" -v minCalls="$minCalls" '
	{
		if(substr($0, 0, 1)!="#") {
			if(length($5)==1 && length($4)==1) {
				refCt=gsub("0/0","0/0",$0)-1; 
				altCt=gsub("1/1","1/1",$0); 
				if(refCt>0 && altCt>0 && (refCt+altCt) > minCalls){
					
					printf $2","$4 
			
					for(i=fsc; i<=NF; i++) {
						split($(i),parts,":")
						GT = parts[1]
						printf "," 
						if(GT=="0/0") {
							printf $4 
						} else if (GT=="1/1") {
							printf $5 
						} 
						else {
							printf "N"
						} 
						
					}
					print "" 
				}
			}
		}		
	}' >> $outFile
}
export -f makeSNPtable

parallel --gnu -j${threads} makeSNPtable $vcf $firstSampleCol $minCalls $outfileBase ::: ${chroms[*]}


	

