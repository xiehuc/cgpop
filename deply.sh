#!/bin/bash

declare -A nprocs
# small set
nprocs['180x120']="360 180  120  96  48  36  24"
nprocs['120x80']="768 672  576  480 384 264 192 96  48"
nprocs['90x60']="1320 1156 990  826 660 552 444 336 252 168 84"
nprocs['60x40']="1182 1065 948  828 708 642 576 468 360 270 180 90"
nprocs['48x32']="1272 1182 984  876 714 552 414 345 276 138 69"
nprocs['36x24']="1304 1245 1186 948 830 712 593 474 356 238 120"

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
	sed -e "s/\$nodes/$(($pe/12))/g" -e "s/\$job/$job/g" -e "s/\$pe/$pe/g" deply.pbs.template > run/$job/deply.$pe.pbs
	cd run/$job
	qsub deply.$pe.pbs
	cd ../..
	done 
done
