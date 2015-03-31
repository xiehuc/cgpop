#!/bin/bash
# a macro script means run all test in single machine, didn't deploy to
# clusters.

declare -A nprocs
# small set                                       #commented would crash
nprocs['180x120']={24..360..3}
nprocs['120x80']={48..768..3}
nprocs['90x60']={84..1320..5}
nprocs['60x40']={96..1188..5}
nprocs['48x32']={72..1272..5}
nprocs['36x24']={96..1320..5}

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
