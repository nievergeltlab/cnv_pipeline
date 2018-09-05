#Step 1: Convert PennCNV/iPattern to .bed format
#Step 2: Add "strand" information to .bed files, where - strand stands for a deletion and + stands for an insertion. 
#This is a trick to allow merge of files for corresponding calls in bedtools, without having to make separate files for duplications and deletions
#Step 3: Merge calls in bedtools
#Step 4:QC subjects
#Step 5: Anneal CNVs
#Step 6: Remove unwanted regions
#Step 7: Rare CNV detection

#Steps 6 and 7 can be done in either order. 

#requires extra input files:
#1)  File listing, one subject per row the ipattern subject name and the filename for the subject that went into PennCNV
#.subjects_list appended to name
#) Conversion sheet consisting of ipattern subject name then GWAS FID and IID
#study_conversion_key.txt
#) family file for GWAS
study='psy4'

cd /home/amaihofer/$study

mkdir temporary_files
 
#####################################
### convert PennCNV results to a BED file (BASH)
#####################################

#PennCNV comes with a script to do this pretty easily

#These files use color code designations to denote which are insertions and deletions

#First concatenate all calls
cat penncnv/calls/*.rawcnv > temporary_files/"$study".rawcnv

#Then use visualize_cnv script to convert to .bed format
#Remove the full paths to files using sed, otherwise it'll screw up the write procedure for R

#Skip the header line with tail -n+2, as we don't want it yet ( it is e.g.: track name="Track: "$study"_PennCNV" description=""$study"_PennCNV" visibility=2 itemRgb="On" )
#Visualizecnv2 is hacked to include nsnps in one of the columns that is otherwise null
/home/amaihofer/PennCNV-master/visualize_cnv2.pl -intype cnv -format bed --track_name "$study"_PennCNV temporary_files/"$study".rawcnv | sed 's/\/\//\t/g' | cut -f 1-3,5- | tail -n+2 > temporary_files/"$study"_allsamples.penncnvbed

#Visualize CNV. pl sucks 
#For QC: additional concatenate all log files
cat penncnv/calls/*.log> temporary_files/"$study".cnvlog

/home/amaihofer/PennCNV-master/filter_cnv.pl  temporary_files/"$study".rawcnv -qclogfile temporary_files/"$study".cnvlog  -qcpassout temporary_files/"$study".qcpass -qcsumout temporary_files/"$study".qcsum -out temporary_files/"$study".goodcnv
#Need to remove those damn paths again 
cat temporary_files/"$study".qcsum | sed 's/\/\//\t/g' | awk '{OFS="\t"}{if(NR==1) print; else print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > temporary_files/"$study".qcsum_fix

#####################################
### Convert iPattern to BED format (BASH)
#####################################

#Convert the iPattern files to the penncnv bed format then just conver them to .bed files
#Thsis format is a little weird because column denotes multiple subjects. Needs to be reparsed
#we'll convert the data into a tab format, using some dummy columns

#iPattern headers: 
#1.	CNV type: gain/loss calls are relative clean calls, cxGain and cxLoss calls are relative noisy regions deemed by the program.
#2.	Chr: chromosome number
#3.	start: the start point of genome coordinate
#4.	end: end end point of genome coordinate
#5.	probe#: how many number of probes are in the region
#6.	on_probe#: how many number of probes used by the algorithms to make CNV call
#7.	clusterIdx: the cluster identification used by the algorithm to make CNV calls; because sometimes the algorithm groups a set of samples together
#8.	gain/loss_score: a suggestive confidence score
#9.	cluster score: a suggestive for the clustered samples
#10.gain_/loss_sample#: how many samples in the cluster
#11.sample_score: how likely the call(s) of the sample, the sample may have two calls nearby
#12.sample_id: which sample this call belongs to
#13.CNV_event_ID
#14. CNVR_ID: identifications used by the algorithm to group the calls..

#bed format chr2	148715037	148758479	CN0003362	100	-	0	0	255,0,0
#1 chr#
#2 start pos
#3 stop pos 
#4 sample name
#5 confidence score
#6 + or -  
#7 null 
#8 null, but now I'm putting in N SNPs, according to col 6
#9 255,0,0 for loss, 0,255,0 for gain

#Currently cxLoss and cxGain are coded as losses and gains
#Currently N snps is counted as on probe # (column 6), rather than columnn 5 (N probes in region). 
#This is important because the filtering step for N probes is done here!
#Concatenate all plate batch files and convert to a .bed format

#Check that resutls directories have produced outputs
ls results_*/

#Concatenate all results files
cat results_*/"$study"_*_all_calls.txt > temporary_files/"$study".ipn

cat temporary_files/"$study".ipn | grep -v \# | awk 'BEGIN{OFS="\t"} {if($1 == "cxLoss" || $1 == "Loss") {copies="-"; color="255,0,0"} else if ($1 == "cxGain" || $1 == "Gain") {copies="+"; color="0,255,0";} print "chr"$2,$3,$4,$12,$8,copies,0,$6,color}'  > temporary_files/"$study"_allsamples.ipnbed

#Tagged results files
cat temporary_files/"$study".ipn | grep -v \# | awk 'BEGIN{OFS="\t"}{print "chr"$2":"$3"-"$4,$0}' > temporary_files/"$study".ipn_tagged

#####################################
### Split BED files by subject and CNV type (R)
#####################################


#Split the files by subject. CNV types (insertion or deletion)
#are noted by the strand column (+/-))
#Insertions and deletions are currently not distinguished by count in this export
#e.g. a 2 deletion CNV will go into the same file as the 1 deletion CNVs

#Split the files by subject and by cnv type
##we split by indels using the color codes
#color codes are:
#1=>'128,0,0', 2=>'255,0,0', 3=>'0,255,0', 5=>'0,255,0', 6=>'0,128,0')
#state 1 = 0 copies
#state 2 = 1 copy 
#state 3 = 2 copies?? (e.g. x chr)
#state 4 = not existing in this data
#state 5 = 3 copies
#state 6 = 4 copies

#load R
R

study='psy4'
library(plyr)

options(stringsAsFactors=F)

tablefun <- function(x, caller)
{
	#Currently spits out into the temporary_files directory
	outfilename <- paste('temporary_files/',x$Sample.ID[1],'.',caller,'bed',sep='')
	write.table(x,outfilename,quote=F,row.names=F,col.names=F,sep='\t')
}

for (caller in c("ipn","penncnv"))
{
 dat <- read.table(paste('temporary_files/',study,'_allsamples.',caller,'bed',sep=''), header=F,sep='\t')

 names(dat)[4] <- "Sample.ID"
 names(dat)[9] <- "DelorDup"
 names(dat)[6] <- "strand"
  
 dat[which(dat$DelorDup %in% c('128,0,0', '255,0,0')),]$strand <- "-"
 dat[which(dat$DelorDup %in% c('0,255,0', '0,128,0')),]$strand <- "+"

 #write a table for each subject
 d_ply(dat, ~ Sample.ID, tablefun,caller=caller)
}

q()
n

#####################################
### Combine CNV calls
#####################################

############################################################
### Intersect and merge CNV calls between programs (BASH)
############################################################

study='psy4'
module load bedtools

nsubs=$(wc -l "$study".subjects_list | awk '{print $1}')

#loop over all subjects, using sample alignment list because iPN and PennCNV have different ways of naming subjects
#Be sure that the subject list file has a linux end of line conversion (ie cant just use excel!)

for i in $(seq 0 1 $nsubs )
 do
  j=$(awk -v linez=$i 'NR==linez{print $1,$2}' "$study".subjects_list)
  #### Note again: we merge on not just overlap, but what we put in the strand column. This will keep opposing calls (e.g. PennCNV insertion, ipAttern deletion) from merging
  ipn_name=$(echo $j | awk '{print $1}')
  penn_name=$(echo $j | awk '{print $2}')
 echo $ipn_name
 echo $penn_name
  #PennCNV and iPattern merge, then sort by Chr and position
  bedtools intersect -s -a temporary_files/"$penn_name".penncnvbed -b temporary_files/"$ipn_name".ipnbed | sort -n -k 1.4,1.5 -k 2,2 > temporary_files/"$ipn_name"_strand.bed_merged_byfile 
  
  #Make same file but with position of what got overlapped from each file and also put in the probe count from ipattern 
  #There is no strong need to have a tagger file. These files going into the merge are at the point where you should have all columns you want to retain!
  bedtools intersect -wa -wb -s  -a temporary_files/"$penn_name".penncnvbed -b temporary_files/"$ipn_name".ipnbed | sort -n -k 1.4,1.5 -k 2,2 | awk 'BEGIN{OFS="\t"}{print $1":"$2"-"$3,$10":"$11"-"$12,$4,$13,$17}'  > temporary_files/"$ipn_name"_strand.overlapbed_merged_byfile 
  
  #modifying probe count to go in here in place of the confidence metric
  #We dont use confience scores so I put the other ID in here
  paste temporary_files/"$ipn_name"_strand.bed_merged_byfile  temporary_files/"$ipn_name"_strand.overlapbed_merged_byfile | awk 'BEGIN{OFS="\t"}{if($7>= $14) probecol=$7; else probecol=$14; print $1,$2,$3,$4,probecol,$6,$7,$14,$9}' > temporary_files/"$ipn_name"_strand.bed_merged_posinfo
  
done

#combine all subjects data into one file

 cat temporary_files/*.bed_merged_posinfo > temporary_files/"$study"_strand.bed_merged

 cat temporary_files/*.overlapbed_merged_byfile > temporary_files/"$study"_strand.overlapbed_merged
 
############################################################
### Remove subjects based on quality control information (R)
############################################################

module load R
R
#Check IQR rules
study='psy4'

library(plyr)

options(stringsAsFactors=F)

#Get chromosome sizes,
#From http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
chr_sizes <- read.table('/home/amaihofer/cnv_reffiles/hg19_chr_sizes.txt', header=F)
names(chr_sizes) <- c("Chrom", "size")

aneuploidy_check <- function(x, chr_sizes,  removepct)
{
	#need a table declaring size of each chromosome in BP 
	#perform the aneuploidy check
	#Namely, for any chromosome,
	#If X% of the chromosome is covered by CNVs, 
	#Subject is aneuploid
	
	#function is meant to be run per subject, per chromosome
	#Get the sum of all CNV lenghts on the chromosome
	for( chromosome in chr_sizes$Chrom)
	{
		x_chr <- subset(x, CHR = chromosome)
		bpdifs <- x_chr$BP1 - x_chr$BP2
		chr_size <- sum(bpdifs)
		if(chr_size/subset(chr_sizes, Chrom == chromosome)$size > removepct)
		{
			return(1)
		}
	}
	return(0)
}


#read the .bed file
infilename <- paste('temporary_files/',study,'_strand.bed_merged',sep='')
dat <- read.table(infilename, header=F,sep='')
names(dat) <- c("CHR","BP1","BP2","SUBJ","LENGTH","STRAND","JUNK1","JUNK2","CNVTYPE")
#Get N CNVs per subject
cnv_count <- ddply(dat, ~ SUBJ,dim)[,c(1,2)]
names(cnv_count)[1] <- "subject"
names(cnv_count)[2] <- "ncnvs"

#Remove subjects with an outlying number of CNV (3 IQR rule)
iqr_rule_cnv <- quantile(cnv_count$ncnvs,.5) + 3*IQR(cnv_count$ncnvs)
remove_count <- cnv_count[cnv_count$ncnvs >= iqr_rule_cnv ,]$subject

#Read in the penncnv generated qc info
infileqc <- paste("temporary_files/",study,'.qcsum_fix',sep="")

qc_file <- read.table(infileqc,header=T)
qc_file$subject <- qc_file$File

#determine what constitutes an outlier in these data and schedule to remove
LRR_iqr <- quantile(qc_file$LRR_SD,.5) + 3*IQR(qc_file$LRR_SD)
baf_iqr <- quantile(qc_file$BAF_SD,.5) + 3*IQR(abs(qc_file$BAF_SD))
waviness_iqr <- quantile(abs(qc_file$WF),.5) + 3*IQR(abs(qc_file$WF))

remove_lrr <- qc_file[qc_file$LRR_SD > LRR_iqr,]$subject
remove_wf <- qc_file[abs(qc_file$WF) > waviness_iqr ,]$subject
remove_baf  <- qc_file[abs(qc_file$BAF_SD) > baf_iqr  ,]$subject


#Perform aneuploidy checks and schedule to remove
aneup_check <- ddply(dat, ~ SUBJ, aneuploidy_check, chr_sizes=chr_sizes, removepct=.1)
names(aneup_check)[2] <- "failed"

remove_aneuploid <- subset(aneup_check, failed == 1)

#remove all failed subjects
kept <- subset(dat ,!(SUBJ %in% remove_count) & !(SUBJ %in% remove_lrr) & !(SUBJ %in% remove_baf) & !(SUBJ %in% remove_wf) & !(SUBJ %in% remove_aneuploid) )
#!(SUBJ %in% remove_dels) &  !(SUBJ %in% remove_dups) & 


outdat <- kept

write.table(outdat, paste('temporary_files/',study,'_sample_qced.merged',sep=''), quote=FALSE, col.names=TRUE, row.names=F,sep="\t")

#write a covariate file of study standardized LRR SD and waviness
covariates <- qc_file
covariates$LRR_SD_standardized <- scale(covariates$LRR_SD,center=TRUE,scale=TRUE)

write.table(covariates, paste(study,"_pcvcovariates.txt",sep=""),quote=F,row.names=F)

############################################################
### Anneal the CNV consensus calls (R)
############################################################

study='psy4'
library(plyr)

cnv_anneal <- function(cnv_matrix, gaplength)
{

#This interrogates every adjacent pair of CNVs to anneal them.
#Input file must be sorted by position!

#The current merging algorithm:
#To merge, the gap between a pair of CNVs must be <= x% of the 
#total distance spanned by the combined length of the CNVs and
#the gap itself.
#Every time the function determines that a merge must be made
#it performs the merge then starts over from the beginning

	i=1
	while (i < dim(cnv_matrix)[1])
	{
		A <- (cnv_matrix[i,]$BP2 - cnv_matrix[i,]$BP1)
		B <- (cnv_matrix[i+1,]$BP1 - cnv_matrix[i,]$BP2)
		C <- (cnv_matrix[i+1,]$BP2 - cnv_matrix[i+1,]$BP1) 

		distance <-  B / (A + B + C)
		
		#if it indeed can be merged, as judging by the gap length,
		#modify the first entry to be the total annealed length, delete the second entry,
		#and then start the loop over from the beginning
  #Also adjust the number of SNPs (length) as the sum of the two segments
		if(distance < gaplength)
		{
			cnv_matrix[i,]$BP2 <- cnv_matrix[i+1,]$BP2
			#cnv_matrix[i,]$LENGTH <- cnv_matrix[i,]$BP1 + cnv_matrix[i+1,]$BP2
   cnv_matrix[i,]$LENGTH <- cnv_matrix[i,]$LENGTH + cnv_matrix[i+1,]$LENGTH
			cnv_matrix <- cnv_matrix[-(i+1),] 
			i=1
		} else i <- i+1
	}
	return(cnv_matrix)
}

dat <- read.table(paste('temporary_files/',study,'_sample_qced.merged',sep=''),header=T,stringsAsFactors=F)


#Note: sometimes one program will call a CNV as being on both chromosomes. Another will call it as being
#only for one chromosome. This makes  for a somewhat discordant annealing: it will. We take the more conservative approach
#of calling these single chromosome CNVs.
cnv_conv <- function(x)
{
	if(x == "128,0,0")
	{
		return("255,0,0")
	}
	if(x == "255,0,0")
	{
		return("255,0,0")
	}
	if(x == "128,0,0,255,0,0")
	{
		return("255,0,0")
	}
	if(x == "0,128,0,0,255,0")
	{
		return("0,255,0")
	}
	if(x == "0,255,0")
	{
		return("0,255,0")
	}
	if(x == "0,128,0")
	{
		return("0,255,0")
	}
}

dat$CNVTYPE <- sapply(dat$CNVTYPE, cnv_conv)

#Notice, this is where the gap length should be set! By default it is 30%
#Function will loop over each subject and CNV type

cnv_fixed <- ddply(dat, ~ SUBJ + CHR + CNVTYPE , cnv_anneal,gaplength=.3)

#revise CNV lengths 
#Commented out. no need for this, also using number of probes as the length metric now
#cnv_fixed$LENGTH <- cnv_fixed$BP2 - cnv_fixed$BP1

write.table(cnv_fixed,paste('temporary_files/',study,'_annealed.bed_merged',sep=''),quote=F,row.names=F,col.names=F,sep="\t")

q()
n


#####################################
### Remove unwanted regions (e.g. centromeres) (BASH)
#####################################

#Notice that segdups and immunoglobulin/tcell stuff are strand specific (-s appended!!)

mkdir calls_merged

#remove centromeres, if a cnv so much as touches the centromere (defined as 0.01% overlap)
 bedtools subtract -a temporary_files/"$study"_annealed.bed_merged -b /home/amaihofer/cnv_reffiles/centromeres.bed -f 0.01 > temporary_files/"$study"_annealed_qc1.bed_merged

#remove telomeres
 bedtools subtract -a temporary_files/"$study"_annealed_qc1.bed_merged -b /home/amaihofer/cnv_reffiles/telomeres.bed > temporary_files/"$study"_annealed_qc2.bed_merged

#remove segdups (50pct overlap)
 bedtools subtract -a temporary_files/"$study"_annealed_qc2.bed_merged -s -b /home/amaihofer/cnv_reffiles/GRCh37GenomicSuperDup_col16.bed -f 0.50 > temporary_files/"$study"_annealed_qc3.bed_merged

#remove immunoglobulin and t cell receptor loci
 bedtools subtract -a temporary_files/"$study"_annealed_qc3.bed_merged -s -b /home/amaihofer/cnv_reffiles/immunoglobulin_tcella_tcellb.bed  -f 0.50 > calls_merged/"$study"_annealed_qcfin.bed_merged

#### 
## To do: add back confidence scores and stuff (probe count) in the dummy columns
 
######################################################
####### convert the file to be run in PLINK (R)
#####################################################
#first we concatenate everything together, then we convert from BED to PLINK ourselves
#note that for color coding, 
R
#state = color code
#1=>'128,0,0', 2=>'255,0,0', 3=>'0,255,0', 5=>'0,255,0', 6=>'0,128,0'
#state 1 = 0 copies 
#state 2 = 1 copy
#state 5 = 3 copies
#state 6 = 4 copies

#set study name
study='psy4'
#run in R
options(stringsAsFactors=F)

#convert color coding to cnv coding
cnv_conv <- function(x)
{
	if(x == "128,0,0")
	{
		return(0)
	}
	if(x == "255,0,0")
	{
		return(1)
	}
	if(x == "128,0,0,255,0,0")
	{
		return(1)
	}
	if(x == "0,128,0,0,255,0")
	{
		return(3)
	}
	if(x == "0,255,0")
	{
		return(3)
	}
	if(x == "0,128,0")
	{
		return(4)
	}
}

dat <- read.table(paste('calls_merged/',study,'_annealed_qcfin.bed_merged',sep=''),header=F)
names(dat) <- c("CHR","BP1","BP2","SUBJ","LENGTH","STRAND","JUNK1","JUNK2","CNVTYPE")
dat$CHR <- gsub("chr","",dat$CHR)
dat$SITES <- dat$LENGTH
dat$TYPE <- sapply(dat$CNVTYPE,cnv_conv)
dat$SCORE <- 0 #Dummt code for confidence score. Not in use right now.

#Need to convert the IDs to a PLINK format using a conversion key

#read in sample IDs conversion sheet (format is: subject genotyping id, subject FID, subject IID)
sample_ids <- read.table(paste(study,'_conversion_key.txt',sep=''), header=T,na.strings=c(NA,"#N/A"))
names(sample_ids) <- c("SUBJ", "FID", "IID")

datA <- merge(dat, sample_ids,by="SUBJ")

dat_exp <- subset(datA, select=c(FID,IID,CHR,BP1,BP2,TYPE,SCORE,SITES))	

write.table(dat_exp, paste("calls_merged/", study,"_sample_qced.cnv",sep=""),quote=F,row.names=F)

q()
n

########################################################
####### convert the file to be run in PLINK (BASH) 
########################################################

#Put all individual studies  into one big file for PLINK
#Multiple FID line repeats because each individual file had a header. This just takes the first instance of the ehader
#cat calls_merged/*.preplink | awk '{if (NR == 1 || $1 != "FID") print}' > calls_merged/GWAS1_2_combined.cnv


#Make a .map file for PLINK
/home/amaihofer/grac/plink-1.07-x86_64/plink --noweb --cnv-list calls_merged/"$study"_sample_qced.cnv --cnv-kb 0 --cnv-make-map --out calls_merged/"$study"_sample_qced

#Make a .fam for plink by using the phenotype datasheet or just copying the one from GWAS
cp pts_"$study"_mix_am-qc.fam  calls_merged/"$study"_sample_qced.fam

#Check for CNVs overlapping eachother (indicating a nonsensical result)
/home/amaihofer/grac/plink-1.07-x86_64/plink --noweb --cfile calls_merged/"$study"_sample_qced --cnv-check-no-overlap --out "$study"_sample_qced_qcmetrics

#Overall burden test
/home/amaihofer/grac/plink-1.07-x86_64/plink --noweb --cfile calls_merged/"$study"_sample_qced --cnv-sites 10 --cnv-kb 20 --cnv-indiv-perm    --mperm 10000 --out temporary_files/burden1test
#breakpoint test
/home/amaihofer/grac/plink-1.07-x86_64/plink --noweb --cfile calls_merged/"$study"_sample_qced  --mperm 10000 --out temporary_files/burden2test
#Gene enrichment test
/home/amaihofer/grac/plink-1.07-x86_64/plink --noweb --cfile  calls_merged/"$study"_sample_qced  --cnv-intersect /home/amaihofer/cnv_reffiles/glist-hg19 --cnv-test-region --mperm 10000 


