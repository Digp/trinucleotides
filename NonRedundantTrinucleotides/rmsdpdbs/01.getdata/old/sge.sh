#$ -S /bin/bash
#$ -cwd
#$ -q "all.q@xesh1"
#$ -N MC
#$ -o std.o 
#$ -e std.e
#$ -pe mpi 16

module purge
#module load R/3.4.4-goolf-1.4.10
module load R/3.6.0-goolf-1.4.10
Rscript get_trinuc.R
