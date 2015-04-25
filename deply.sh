#!/bin/bash

declare -A nprocs
# small set                               
nprocs['180x120']="24 36 48 96 120 180 360"
nprocs['120x80']="96 192 264 384 480 576 672 768"
nprocs['90x60']="168 252 336 444 552 660 8282 996 1164 1248 1320"
nprocs['60x40']="96 180 276 360 468 576 648 708 828 948 1068 1188"
nprocs['48x32']="72 144 348 420 552 720 876 984 1188 1272"
nprocs['36x24']="96 180 240 360 480 600 720 840 960 1200 1260 1320"

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
