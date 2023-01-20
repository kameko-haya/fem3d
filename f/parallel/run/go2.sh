#@$-r go2
#@$-q lecture
#@$-N 1
#@$-J T8
#@$-e err
#@$-o b08.lst
#@$-lM 27GB
#@$-lE 00:15:00
#@$-s /bin/sh
#@$

cd $PBS_O_WORKDIR
mpirun ./n2.sh ./sol

exit
