cgpop
=====

a mirror of cgpop https://code.google.com/p/cgpop/

use wget to download other data file (too big to put into github):

   wget https://cgpop.googlecode.com/svn/trunk/data/topography.0.1 -O data/topography.0.1

download http://www.cs.colostate.edu/hpc/cgpop/cgpoptiles.tgz to extract to
data dir, or you can use cginit to generate these files according to pdf file
in doc

provide new `deply.sh` to run cgpop. it would write pbs under run folder. you
must run at `run` folder to run cgpop correctly
