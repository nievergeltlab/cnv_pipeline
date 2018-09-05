study=psy4

wd=/home/amaihofer/"$study"/penncnv
mkdir $wd
cd $wd

mkdir calls

intensity_dir=/oasis/scratch/comet/amaihofer/temp_project/"$study"/intensities/
intensitypcv_dir=/oasis/scratch/comet/amaihofer/temp_project/"$study"/intensities_penncnv/

mkdir $intensitypcv_dir

#Do once: Get GC content info for hg19 (if this is one of the cogend studies, hg18)
 # wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gc5Base.txt.gz
 # sort -k2,2 -k3,3n ../gc5Base.txt > gc5Base.txt.sorted

#Note: PennCNV installed in: /home/amaihofer/PennCNV-master

##Make PFB file


#Split data by columns (my way)
#MAKE SURE columns match, peolpe have a habit of fliping SNP and sample name!!
#Note how I put to just remove Nan...
for files in $(ls $intensity_dir )
do
 awk -v datname=$files 'NR>10{print $1,$9,$10}' "$intensity_dir"/$files | sed 's/\r//g' | awk 'BEGIN{OFS="\t"} {if(NR==1) print "Name","datname.B Allele Freq", "datname.Log R Ratio" ; else print $1,$2,$3}' > "$intensitypcv_dir"/$files
done


#Get list of SNP positions
#awk -F"\t" 'BEGIN{OFS="\t";} NR>9{print  $2,$4,$5}' "$intensity_dir"PTSD_9900 > /home/amaihofer/"$study"/penncnv/"$study"snppostions.txt #If no probe file
#if probe file just use that
echo -e "SNP Name\tChr\tPosition" > snpposheader.txt
cat snpposheader.txt  <(awk 'BEGIN{OFS="\t"}NR>11{print $1,$5,$6}' "$intensity_dir"/$files ) > "$study"snppostions.txt 

#List random set of intensity files for generating PFB file
ls "$intensitypcv_dir"* | sort -R | head -n 300 > /home/amaihofer/"$study"/penncnv/"$study"pfb_input.txt

#Make PFB file 
/home/amaihofer/usr/perl-5.14.2/bin/perl /home/amaihofer/PennCNV-master/compile_pfb.pl --listfile  "$study"pfb_input.txt --snpposfile   /home/amaihofer/"$study"/penncnv/"$study"snppostions.txt --out "$study".pfb

#Make GC model file
/home/amaihofer/PennCNV-master/cal_gc_snp.pl gc5Base.txt.sorted /home/amaihofer/"$study"/penncnv/"$study"snppostions.txt --out "$study"_gcmodel

#CNV detection test script
for files in $(ls "$intensitypcv_dir" | tail -n3 ) # 202 204 206 446 668 need redone
do
 # REPLACE PTSD_9900 with 'files' to do all files, otherwise this is just a testrun
 /home/amaihofer/usr/perl-5.14.2/bin/perl /home/amaihofer/PennCNV-master/detect_cnv.pl -test -hmm /home/amaihofer/PennCNV-master/lib/hhall.hmm \
 -pfb "$study".pfb "$intensitypcv_dir"$files -conf -gcmodel  "$study"_gcmodel -log calls/"$files".log \
 -out calls/"$files".rawcnv
done

#List all intensities files
 ls $intensitypcv_dir | grep -v sample_ids > "$study"_samples.txt
 
  
#Split list of samples into sub lists, to be sent to job script
 split -d  -l 250 /home/amaihofer/"$study"/penncnv/"$study"_samples.txt "$study"split
 
#Fix the naming of files to remove padding, so the job array works. It's a little bit fucked up right now
for x in $(ls "$study"split*) 
do
mv -i -- "$x" "${x%/*}.$(printf %05d $((10#${x##*.}+1)))";
done
#Give correct N of arrays..

sbatch --array=1-3 ../call_cnv_serial.qsub -e /home/amaihofer/"$study"/penncnv/errorlog1




#Experimental stuff with using more than 1 proc.
 # nodeuse=24
 # sbatch --array=1 call_cnv.qsub 
 # sbatch call_cnv_serial.qsub 


