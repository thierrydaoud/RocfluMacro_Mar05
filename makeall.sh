module purge
module load gcc/12.2.0 openmpi/4.1.5
module list
cd libpicl
make clean ; make 
cd .. 
make clean ; make RFLU=1 PICL=1 FOLDER=1 -j16
ls --color=auto
