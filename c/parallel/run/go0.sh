#@$-r go
#@$-q debug
#@$-N 1
#@$-J T8
#@$-e err
#@$-o t01.lst
#@$-lM 27GB
#@$-lE 00:05:00
#@$-s /bin/sh
#@$


cd $PBS_O_WORKDIR
mpirun ./sol

exit
