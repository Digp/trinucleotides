#!/bin/bash

for i in $(ls | grep _out); do 


sbatch <<EOF 
#!/bin/tcsh 
#SBATCH -J $i
#SBATCH -c 1
#SBATCH -p MPI

cd /orozco/projects/PDB.datamining/rmsd_data3/$i; 
wc -l *pdb > /orozco/homes/pluto/ematamoros/Desktop/Caca/$i.txt

EOF

done

