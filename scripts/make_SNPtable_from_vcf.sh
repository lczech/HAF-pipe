#!/bin/bash


# ==================================================================================================
#      Usage
# ==================================================================================================

usage()
{
    echo "usage: make_SNPtable_from_vcf.sh [-v --vcf] [-c --chrom]
    optional:
    [ -s --snptable (default: $(echo $vcf | sed 's/.gz$//' | sed 's/.vcf//').$chrom.snptable)]
    [ -f --firstSampleCol (default 10) ]
    [ -m --mincalls (default 2) ]
    [ -u --subsetlist (default none)]
    [ -k --keephets ]
    [ -t --threads (default 1) ]
    [ -h --help ]
    ** note that the vcf may be gzipped
    ** only biallelic sites with at least one ref and one alt call will be recorded
    ** use --mincalls to limit to sites with a higher number of called alleles
    ** subset list is a 1-column list of names of founder haplotypes to keep in the snptable;
       if not supplied, all founders will be included
"
}
if [ $# -lt 1 ]; then
    usage
    exit 1
fi

# ==================================================================================================
#      Command Line Arguments
# ==================================================================================================

firstSampleCol=10
threads=1
keephets=0
mincalls=2
subsetlist=none

while [ "$1" != "" ]; do
    case $1 in
        -v | --vcf )            shift
                                vcf=$1
                                ;;
        -c | --chrom )          shift
                                chrom=$1
                                ;;
        -f | --firstSampleCol ) shift
                                firstSampleCol=$1
                                ;;
        -m | --mincalls )       shift
                                mincalls=$1
                                ;;
        -k | --keephets )
                                keephets=1
                                ;;
        -u | --subsetlist )     shift
                                subsetlist=$1
                                ;;
        -s | --snptable )       shift
                                snptable=$1
                                ;;
        -t | --threads )        shift
                                threads=$1
                                ;;
        -h | --help )           usage
                                exit 1
                                ;;
        * )                     echo "unknown flag $1"
                                usage
                                exit 1
    esac
    shift
done

if [ -z $snptable ]; then
    snptable=$(echo $vcf | sed 's/.gz$//' | sed 's/.vcf//').$chrom.snptable
fi

# ==================================================================================================
#      Main
# ==================================================================================================

echo "making snptable for chrom $chrom from $vcf, starting from column $firstSampleCol"
maindir=$(dirname $0)/..

# In the following, we use gunzip instead of zcat, for portability with MacOS,
# see https://serverfault.com/a/704521 for details.

##PRINT HEADER
gunzip -c $vcf | head -1000  | grep "#CHROM" | head -1 | \
awk -v fsc=$firstSampleCol '{for(i=fsc;i<=NF;i++){printf $i","}}' | \
sed "s/^/${chrom},Ref,/" | sed 's/,$/\n/' > $snptable

gunzip -c $vcf | grep -P '^'$chrom'\t' | \
awk -v fsc="$firstSampleCol" -v mincalls="$mincalls" -v keephets="$keephets" '
BEGIN{
    baseCodes["AG"]="R";baseCodes["GA"]="R"
    baseCodes["CT"]="Y";baseCodes["TC"]="Y"
    baseCodes["CG"]="S";baseCodes["GC"]="S"
    baseCodes["AT"]="W";baseCodes["TA"]="W"
    baseCodes["GT"]="K";baseCodes["TG"]="K"
    baseCodes["AC"]="M";baseCodes["CA"]="M"
}
{
    if( length($4)==1 && length($5)==1 && ( $7 == "PASS" || $7 == "." )) {

        gt_field=0
        split($(fsc-1),field_ids,":")
        for(ii=1;ii<=length(field_ids);ii++){
            if(field_ids[ii]=="GT"){gt_field=ii}
        }
            if(gt_field==0){
                print "ERROR: Genotype (GT) not found in VCF info fields! Exiting."; exit
            }

            genotypes=""
            for(ii=fsc;ii<=NF;ii++){
                split($ii,fields,":");
                gsub("\\|","/",fields[gt_field]);
                genotypes=genotypes","fields[gt_field]
            };

        refCt=gsub("0/0","0/0",genotypes);
        altCt=gsub("1/1","1/1",genotypes);
        hetCt=gsub("0/1","0/1",genotypes);
        missingCt=gsub("\\./\\.","./.",genotypes);

        if((hetCt>0 || (refCt>0 && altCt>0)) && ((refCt+altCt+hetCt) >= mincalls)){

            printf $2","$4
            split(genotypes,gt_array,",")
            for(ii=2; ii<=length(gt_array); ii++) {

                GT = gt_array[ii]
                printf ","
                if(GT=="0/0") {
                    printf $4
                } else if (GT=="1/1") {
                    printf $5
                } else if (keephets>0 && ((GT=="0/1") || (GT=="1/0"))) {
                    printf baseCodes[$4$5]
                }
                else {
                    printf "N"
                }

            }
            print ""
        }
    }
}' >> $snptable

if [ ! "$subsetlist" == "none" ]; then
    ${maindir}/scripts/subset_SNPtable.sh \
    -s $snptable \
    -o ${snptable}.subset \
    -f $subsetlist \
    -m $mincalls ;
    mv ${snptable}.subset $snptable
fi

${maindir}/scripts/count_SNPtable.sh $snptable
${maindir}/scripts/prepare_SNPtable_for_HAFcalc.sh $snptable

echo "SNP table written to: $snptable"
