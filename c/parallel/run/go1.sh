#@$-r go1
#@$-q debug
#@$-N 1
#@$-J T4
#@$-e err
#@$-o t08c.lst
#@$-lM 27GB
#@$-lE 00:05:00
#@$-s /bin/sh
#@$

cd $PBS_O_WORKDIR
mpirun ./n1.sh ./sol

exit
