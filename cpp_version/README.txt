//on compute Canada
module load CCEnv nixpkgs/16.09 gcc/7.3.0 openblas fftw root/6.08.02
export ROOTSYS=/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/root/6.08.02
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export PATH=$ROOTSYS/bin:$PATH

make

run with ./reco to view graphs in real time in ROOT Application

OR

run with ./reco 0 to disable ROOT Application




C++14 is used with ROOT 6.08.02

