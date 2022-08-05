#!/bin/bash
#SBATCH --job-name=OM_delete_files                														
#SBATCH --cpus-per-task=1                  													
#SBATCH --mem=1G              														
#SBATCH --account=penny
#SBATCH --qos=1week          																
#SBATCH --error=/scicore/home/penny/GROUP/M3TPP/<your_exp_name_here>/errors0.txt
#SBATCH --out=/scicore/home/penny/GROUP/M3TPP/<your_exp_name_here>/out0.txt

#############################

# set up
ml purge
cd /scicore/home/penny/GROUP/M3TPP/<your_exp_name_here>
mkdir empty/

# group folders to delete
ls | grep 'err\|base\|scenario' >> del.txt

# delete files and folders
while read filename; do
  rsync -a --delete empty/ $filename/
  echo "Deleting $filename"
  rmdir $filename/
done < del.txt

rm param_tab_*.txt

# clean up
rm del.txt
rmdir empty/
echo "Clean up complete"