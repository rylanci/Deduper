#!/usr/bin/env bash 

#SBATCH --account=bgmp          ### SLURM account which will be charged for the job
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=dedup          ### Job Name
#SBATCH --output=dedup.out         ### File in which to store job output
#SBATCH --error=dedup.err          ### File in which to store job error messages
#SBATCH --time=0-01:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job (usually 1)
#SBATCH --ntasks-per-node=1     ### Nuber of tasks to be launched per Node (usually 1)
#SBATCH --cpus-per-task=1       ### Number of cpus (cores) per task


# activate our python 3 env
conda deactivate  
conda activate bgmp_py3

# Load modules
module load samtools


#Run samtools view and sort on our input sam file
#Output sorted version
samtools view -S -b /projects/bgmp/shared/deduper/Dataset3.sam > temp.bam
samtools sort temp.bam -O sam -o sorted_Dataset3.sam

#Run our deduplicating python script
/usr/bin/time -v ./Deduper_Code.py -i sorted_Dataset3.sam -j Dataset3_deduped.sam -u STL96.txt   




