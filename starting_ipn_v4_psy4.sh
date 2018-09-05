#Extract Ipattern
 tar xvzf pbs_ipn_0.581.tar.gz

  #Notes: for you only group is ddp195, nodes = 1 , proc_per_node = 1, compute partition , for  1 hour.  But this is risky.  
 #And if you run for an hour you will be charged for the full node 24 SUS,  even if you run a serial job and only run on one core.  
 #This means I can make it use less procs to save space.. but check cpu use first, its possible that the apps running ARE multicore!
 
#Setup 1) You have to make some changes to the iPattern scripts in order for it to be compatible with Comet:

#lines 37-40 of ipn/run_ipattern.py: Where there are '$', replace with 'PBS'
 #Reason: We explicitly want to say it's the PBS scheduling system. This will allow comet to conver the parameters to SLURM paramters
 
 
#Line 40 of ipn/run_ipattern.py: It says:   print >> outf, "#PBS -l h_vmem=8G"    . Replace it with this:  print >> outf, "#PBS -l vmem=8G"
 #Reason: vmem command is compatible with comet, but h_vmem is not
 
#Line 40 of ipn/run_ipattern.py: set walltime to #SBATCH -t 00:45:00, 
#reasonthe program runs quicker than expected 

 
 #Line 41 of ipn/run_ipattern.py: Add     print >> outf, "#PBS -l ppn=1" 
 #Reason: Save CPU usage by setting that only 1 proc is used
 
 #line 180 of run_ipattern.py: comment out the sys.exit(1)
 #Reason: itll quit even if nothing is wrong, due to the fact that chr 21 p arm is often not genotyped!!
 
#Line 275 of ipn/run_ipattern.py: One the lines before it says 'copy ( os.path.join(ipn_path, 'known.cnvr.txt')...', add in these two lines of code:
 
#if os.path.exists(os.path.join(tmpdir, 'known.cnvr.txt')):
# os.remove(os.path.join(tmpdir, 'known.cnvr.txt'))

 #Reason: For some reason when this file link gets created, it only has permission to read. If it tries to create the link again, it can't overwrite. and will crash. This deletes the link if it exists already, so it can create it again

#Line 22 of ipnlib/ipn_pbs_qsub.py: Change sleep time from 30 to 200
#Line 452 of ipnlib/ipn_pbs_qsub.py: Change sleep time from 60 to 120
#Line 714 of ipnlib/ipn_pbs_qsub.py: Change sleep time from 60 to 120
 #Reason: It polls the job status system too quickly, making it appear as if some jobs hadn't been submitted, which causes iPattern to fail.
 
#line 262 of ipnlib/ipn_pbs_qsub.py: replace 24 hours with 2 hours or less, e.g. print >> f, '#PBS -l walltime=0:45:00' 
 #Reason: set walltime down from 24 hours to 45. 24 is absurdly long and will delay you in the queue. If you run out of time, just up this value.
#line 263 of ipnlib/ipn_pbs_qsub.py:     print >> f, '#SBATCH --nodes=1' and   print >> f, '#SBATCH --ntasks-per-node=1'
 #Reason: using full node otherwise, cumbersome!!
        
    
 
 #Line 713-714 of ipnlib/ipn_pbs_qsub.py:	says "
 #if not has_jobs_running (waiting_files.keys(), current_jobs) and \
	#					cur >= len(job_files): break
 
 #Line 61 of ipn/iPattern.Runner.R , comment.char = ""
 #Reason: If you have # characters in SNP names, which happens sometimes, it'll crash. This makes it so the hashes are read as characters and not comments.

##Line 94 of common/split_by_gender.py : Comment out sys.exit(1). 
#Reason: The program will quit if it cant split by gender
#Setup2) Open .bashrc file in your home directory. Export iPattern paths in .bashrc file, e.g. 

# export IPNBASE="/home/amaihofer/ipn_0.581"
# export PATH=$PATH:"$IPNBASE/ipn"
# export PATH=$PATH:"$IPNBASE/preprocess/ilmn"
# export PATH=$PATH:"$IPNBASE/preprocess/affy"
# export PYTHONPATH=$PYTHONPATH:"$IPNBASE/ipnlib"

#But without the hash tags included. Replace amaihofer with your username.


#Setup 3) Use shortcuts to tell iPattern that you want to use pbs job scripts , and NOT SGE scripts.

 cd "$IPNBASE"/ipnlib 
 ln -s ipn_pbs_qsub.py ipn_qsub.py #you may have to delete existing ipn_qsub.py link for this to work

#Setup 4) Install ppc R library (only has to be done one time)

 module load R
 wget http://statweb.stanford.edu/~tibs/PPC/Rdist/ppc_1.02.tar.gz
 R
 install.packages('ppc_1.02.tar.gz',repos=NULL)


#iPattern test run code example:

module load R #important to load R every time prior to running ipattern, otherwise it won't find R and fail.
module load python #important to load python (ipattern uses numpy) prior to running ipattern, otherwise it will fail because it can't load numpy.

#Note: For running the job, the scratch space should be /oasis/scratch/comet/$USER/ , as that is the designated scratch space

#There are two ways of running iPattern:
#For whichever way you run it, a bunch of shortcut files will be made in whatever directory you launch from, so I like to start within a specific folder, e.g.
study=psy4
 mkdir /home/amaihofer/$study
 cd  /home/amaihofer/$study
 
 
#Split data by subject
#make_plates.pbs lets you do this, only necessary if data aggregated


#List all samples on a given plate (only necessary if data divided by plates
# for plate in  {1..8}
# do
 # echo "/oasis/scratch/comet/amaihofer/temp_project/safr/intensities/GS_Plate"$plate"_FinalReport.txt" > plate"$plate"_tobesplit.txt
 # awk 'NR>=11{print $2}' /oasis/scratch/comet/amaihofer/temp_project/safr/intensities/GS_Plate"$plate"_FinalReport.txt | sort -u > plate"$plate"_samples.txt
 # #paste plate"$plate"_samples.txt plate"$plate"_samples.txt > plate"$plate"_samples.txt.rename
# done


#Put intensity files in the temp folder 
#Use the script to make intensity file lists. Will attempt to make plates with 200 people
#Rscript make_intensitylists_v2.r "$study"_cnv.csv  $study


#Method 1: From the shell, as a parameterized command.
#Here I set variables to paths
#Note: Do not have carriage returns (windows style returns) in this file, the bad sample file, or data list file. run dos2unix on them first if they came from excel
 gender_file=/home/amaihofer/"$study"/gender_file.txt #Path to gender file.
 bad_sample_file=/home/amaihofer/"$study"/bad_samples.txt #Path to bad sample file
 data_file_list=/home/amaihofer/"$study"/"$study" #List of intesity data files. Should be named for doing a loop!
 split=no #Split data? set to no if this is not needed
 probe_file=xxxx # #Path to probe file. Set to xxxx if unused /home/amaihofer/safr/safr_probefile.txt i dont need it, this program is also fucked trying to use a separate one!
 batch_file=xxxx #path to batch file. Set to xxxx if unused
 output_dir=/home/amaihofer/"$study"/results #Path to output results 
 experiment="$study"  #experiment signature
 tempdir=/oasis/scratch/comet/$USER/temp_project

 if [ $split != "no" ]
 then
  split_command="--split"
 fi

 if [ $batch_file != "xxxx" ]
 then
  batch_command="--b $batch_file"
 fi

 for i in {2..4} # {11..39} # $(ls *intensities* | wc -l | sort) # for each intensity file set, run this command..
 do
  ilmn_run.py -g $gender_file -m $bad_sample_file -f "$data_file_list"_"$i" -x $experiment"_$i" $batch_command $split_command --temp-prefix-directory-name $tempdir --dest-prefix-directory-name $tempdir --call-prefix-directory-name $tempdir  --out "$output_dir"_"$i"
 done
 
#I didn't make a conf file for this project yet
#Method 2) To use configuration files instead of supplying paramters to the command line
##ilmn.sh /home/amaihofer/test.conf
