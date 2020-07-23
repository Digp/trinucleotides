## tleap
#rootname=$(basename $1 .pdb)
rootname=$1
echo $rootname
echo "source leaprc.DNA.bsc1" > tleap.in
echo "source leaprc.RNA.OL3" >> tleap.in
echo "x=loadpdb $1.pdb" >> tleap.in
echo "savepdb x $rootname""_leap.pdb" >> tleap.in
echo "quit" >> tleap.in

tleap -f tleap.in

rm tleap.in
