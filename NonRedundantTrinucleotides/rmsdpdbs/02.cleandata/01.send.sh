#!/bin/bash

#for dir in $(ls -d /orozco/projects/PDB.datamining/rmsd_data3/* | grep /[ACGU]-[AGCU]-[ACGU]$); do
dir=$1

sbatch <<EOF
#!/bin/bash
##### ntasks=NUMBER OF CORES>=NUMBER OF GPUS
#SBATCH -c 1
##### select a partition MPI
#SBATCH -p MPI
###### job name
#SBATCH -J $dir.2
##### error file
#SBATCH -e err.$dir.2.e
###### output file
#SBATCH -o out.$dir.2.o

hostname

module purge
module load R/3.4.4-foss-2018a-X11-20180131

Rscript remove_h_p2.0.R -p $dir
EOF

#done
