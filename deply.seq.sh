#!/bin/bash
# a macro script means run all test in single machine, didn't deploy to
# clusters.

declare -A nprocs
# small set                                       #commented would crash
nprocs['180x120']=$(seq 24 3 360)
nprocs['120x80']=$(seq 48 3 768)
nprocs['90x60']=$(seq 84 5 1320)
nprocs['60x40']=$(seq 96 5 1188)
nprocs['48x32']=$(seq 72 5 1272)
nprocs['36x24']=$(seq 96 5 1320)

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
