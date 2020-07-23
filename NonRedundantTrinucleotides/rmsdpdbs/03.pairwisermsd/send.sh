#!/bin/bash

for i in $(ls ../../../../rmsd_data3/ | grep _out$); do

a=$(basename $i _out)
echo $a

##SBATCH -p MPI
##sbatch --exclude=h-13,h-14,h-39,h-40 <<EOF
sbatch <<EOF
#!/bin/tcsh
#SBATCH -J $a
#SBATCH -c 1
#SBATCH -p MPI

module purge
module load R/3.4.4-foss-2018a-X11-20180131

Rscript calculatermsd.R --path '../../../../rmsd_data3/'$i'/' --file1 PDB_ANALYSIS.txt --file2 PDB_VAL.txt --ofile $a'.txt'
EOF

done
#
