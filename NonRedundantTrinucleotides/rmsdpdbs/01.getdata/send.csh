#!/bin/tcsh

sbatch <<EOF
#!/bin/tcsh
#SBATCH -J trinuc
#SBATCH -c 32
#SBATCH -p XESH

module purge
module load R/3.6.0-goolf-1.4.10
Rscript get_trinuc_v2.R

end 

