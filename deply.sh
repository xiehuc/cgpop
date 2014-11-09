#!/bin/bash

declare -A nprocs
# small set
nprocs['180x120']="360 180 120 96 48 36 24"
nprocs['120x80']="768 384 264 192 96"
nprocs['90x60']="1320 660 444 336 168"
nprocs['60x40']="1416 948 708 576 360"
nprocs['48x32']="1452 1092 876 552"
nprocs['36x24']="1896 948"

#huge set
#nprocs['180x120']="360 180 120 96 48 36 24"
#nprocs['120x80']="768 384 264 192 96"
#nprocs['90x60']="1320 660 444 336 168"
#nprocs['60x40']="2832 1416 948 708 576 360"
#nprocs['48x32']="2172 1452 1092 876 552"
#nprocs['36x24']="2520 1896 948"
for job in ${!nprocs[*]}
do
   for pe in ${nprocs[$job]}
	do
   echo "$job.$pe"
	#sed -e "s/\$nodes/$(($pe/12))/g" -e "s/\$job/$job/g" -e "s/\$pe/$pe/g" deply.pbs.template > run/$job/deply.$pe.pbs
	#cd run/$job
	#qsub deply.$pe.pbs
	#cd ../..
	done 
done
