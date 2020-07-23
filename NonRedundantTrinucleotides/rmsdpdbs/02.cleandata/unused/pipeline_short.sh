#! /bin/bash

#for dir in $(ls | grep '-' | grep -v leap); do

    dir=$1
    pwd
    mkdir $dir'_leap'
    tmp=$dir'_tmp'
    mkdir $tmp
    cd $tmp
    ln -s ../tleap.sh 
    ln -s ../remove_h_p.R 
    
    for i in $(ls ../$dir);
    do
    	echo $i
    	out=$(basename $i .pdb)
    
    	ln -s '../'$dir'/'$out'.pdb' .
    	./remove_h_p.R -f $out'.pdb'

    	./tleap.sh $out
    	mv $out'_leap.pdb' '../'$dir'_leap/'$i
        rm $out'.pdb'    
    	#cd ..
    	#rm -rf $out
    
    done
    cd ..
    rm $tmp
#done
