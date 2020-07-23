#!/bin/bash

#for i in $(ls *txt | grep '-'); do
for i in $(ls *txt | grep '-'); do
    a=$(basename $i .txt)
    echo $a
    #SBATCH -p MPI
    #sbatch --exclude=h-13,h-14,h-39,h-40 <<EOF
sbatch <<EOF
#!/bin/tcsh
#SBATCH -J $a
#SBATCH -c 1
#SBATCH -p MPI

module purge
module load R/3.4.4-foss-2018a-X11-20180131

Rscript clustering.R -i $i -o $a'representatives.txt'
EOF

done


