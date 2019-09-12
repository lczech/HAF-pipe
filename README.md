# HAF-pipe 🏂

HAF-pipe is a bash- and R-based pipeline to calculate haplotype-inferred allele frequencies from pool-seq data and founder SNP profiles.

## Paper



> Tilk et al. High accuracy haplotype-derived allele frequencies from ultra-low coverage pool-seq samples, bioRxiv, 2018. https://www.biorxiv.org/content/early/2018/01/11/244004  

## Dependencies



1. HARP (download from https://bitbucket.org/dkessner/harp )
2. R version >= 3.2 (in your path) + library: data.table 
3. tabix and bgzip (http://www.htslib.org/download/)

## Tasks
*	1 - make SNP table from VCF
*	2 - impute SNP table
*	3 - infer haplotype frequencies
*	4 - calculate allele frequencies



## Usage



```
HAFpipe_wrapper.sh	[ -t --tasks ]		comma-separated list of tasks to run 
			[ -l --logfile ]    	name of file to write log to
						   - default: HAFpipe-log.`date +%Y-%m-%d.%H%M%S` 
			[ -d --scriptdir ]  	directory in which HAF-pipe scripts are located; tasks:1,2,3,4
						   - default: directory of wrapper script 
			[ -o --outdir ]     	output directory; tasks:1,3,4 
						   - default: current directory
			[ -v --vcf ]        	a multi-sample VCF file of each individual founder in the starting population, to be converted to snp table; tasks:1
			[ -c --chrom ]      	name of chromosome to extract from vcf; tasks:1
			[ -s --snptable ]   	snp table to use for calculating haplotype and allele frequencies; tasks:2,3,4 
						   - columns are: position, ref allele, followed by genotypes of the each individual founder 
						   - See 'simulations/99.clean.SNP.HARP.segregating.gz' included in this repository for examples 
                            			   - will be overwritten if task 1 is run in conjunction with other tasks 
			[ -m --method ]     	method to use for imputation in task 2 or file extension for task 3; tasks:2,3
        		    			   - for task 2, method must be one of: 
                            			      - 'simpute' (simple imputation) 
                            			      - 'npute' (see Roberts et al., 2007 - doi:10.1093/bioinformatics/btm220 ) 
                            			   - for task 3, method can be: 
                            		              - 'none' or not specified (default; the original [snptable] with potential missing calls will be used to infer haplotype frequencies) 
                            			      - any other string (will look for the file [snptable].[method] to infer haplotype frequencies) 
			[ -n --nsites ]     	number of neighboring sites to use with 'npute' imputation; tasks:2
                            			   - default: 20 
			[ -b --bamfile ]    	name of bamfile generated by pooled sequencing of sampled chromosomes from an evolved population; tasks:3,4 
			[ -r --refseq ]     	reference sequence in FASTA format; tasks:3
			[ -e --encoding ]   	base quality encoding in bam files; tasks:3
                            			   - 'illumina' (Phred+64; default) 
                            			   - 'sanger' (Phred+33) 
			[ -g --generations ] 	number of generations of recombination; used to calculate window size for haplotype inference; tasks=3
			[ -a --recombrate ] 	recombination rate used to calculate window size for haplotype inference; tasks=3
						   - default: .0000000239 (avg D. melanogaster recomb rate) 
			[ -q --quantile ]   	quantile of expected unrecombined segment distribution to use for determining haplotype inference window size; tasks:3
						   - default: 18 
			[ -w --winsize ]    	user-defined window size (in kb) for haplotype inference; tasks:3
                            			#(overrides -g and -a) 
			[ -h help ]         	show this help screen
```

## Example

Creates a SNP table from a VCF of founder genotypes and uses it to estimate HAFs from a pooled sequencing sample.  
In this example, inferred haplotype frequencies will be written to ```exptData/sampledInds_cage2_gen15.mapped_dm5.freqs```
and HAFs will be written to ```exptData/sampledInds_cage2_gen15.mapped_dm5.afSite```  <br>

``` 
HAFpipe_wrapper.sh -t 1,2,3,4 \
-v refData/founderGenotypes.vcf \
-c 2L \
-s refData/founderGenotypes.segregating.snpTable \  ## will be created
-m simpute \
-b exptData/sampledInds_cage2_gen15.mapped_r5.39.bam \
-e sanger \
-g 15 \
-r refData/dmel_ref_r5.39.fa \
-o exptData 

```

## Simulations

To simulate recombination and pooled sequencing data from a panel of sequenced founders, see the scripts ```HAFpipe-sim.run_forqs.sh``` and ```HAFpipe-sim.simulate_poolseq.sh``` in the simulations folder. 
* You can use the first script to simulate recombination for an experimental population of any size initiated from any sequenced founder set for any number of generations with any number of selected sites with any selection coefficient.  
* You can use the second script to 'sample' any number of individuals from the recombined population at any generation and simulate pooled sequencing (150bp paired-end reads) of these individuals at any coverage with any sequencing error rate.  This script will also record the true allele frequencies of all segregating sites in the set of sampled individuals. Compare HAFs calculated from the simulated sequencing data to these true allele frequencies to assess effective coverage.  

(Note that more generations, larger pool sizes, and higher coverage values will all increase run time)
