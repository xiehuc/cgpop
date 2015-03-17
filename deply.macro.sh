#!/bin/bash
# a macro script means run all test in single machine, didn't deploy to
# clusters.

declare -A nprocs
# small set                                       #commented would crash
nprocs['180x120']="360 180  120  96  48  36  24"
nprocs['120x80']="768 672  576  480 384 264 192 96" #48
nprocs['90x60']="1320 1248 1164 996  828 660 552 444 336 252 168" #84
nprocs['60x40']="1188 1068 948  828 708 648 576 468 360 276 180 96"
nprocs['48x32']="1272 1188 984  876 720 552 420 348 144 72" #276
nprocs['36x24']="1320 1260 1200 960 840 720 600 480 360 240 180 96" #120

#huge set
#nprocs['180x120']="360 180 120 96 48 36 24"
#nprocs['120x80']="768 384 264 192 96"
#nprocs['90x60']="1320 660 444 336 168"
#nprocs['60x40']="2832 1416 948 708 576 360"
#nprocs['48x32']="2172 1452 1092 876 552"
#nprocs['36x24']="2520 1896 948"

mdlist="mpi2s1D"

if [ -z "$ARCHDIR" ]; then
   echo 'please set $ARCHDIR Environment'
   exit 0
fi

for job in ${!nprocs[*]}
do
   for pe in ${nprocs[$job]}
   do
      cd run/$job

      for md in $mdlist
      do
         for cgpop in cgpop cgpop_db
         do
            CGPOP=../../$md/$cgpop.$ARCHDIR.$job
            if [ -e $CGPOP ];then
               echo "running $CGPOP with $pe process" >&2
               # /usr/bin/time, not bash time keyword. -p is a familar output format
               PROFILING_OUTDIR="$cgpop.$ARCHDIR.$pe" MPI_RANK=0 MPI_SIZE=$pe time -p $CGPOP
            fi
         done
      done

      cd ../..
   done 
done
