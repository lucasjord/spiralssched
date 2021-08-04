#!/bin/bash

cd /home/observer/correlations2/spiralssched/ && git pull 

export class=$(echo $1 | sed 's/.$//')

echo "Newsmerd for slogit"
rsync -Puv /home/observer/correlations2/spiralssched/vex/$1.vex observer@newsmerd:/vlbobs/lba/$class/$1.vex


echo "Hobart12"
# Copy vex to pcfshb home area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfshb:~/
# drudg the vex file there
ssh oper@pcfshb 'cd /usr2/sched; echo -e "hb\n3\n0\n"|/usr2/fs/bin/drudg "$1.vex'
# Run the snp2flex script
ssh oper@pcfshb snp2flex.py $1hb.snp

# Copy vex file over the pcfs-2hb sched area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2hb:/usr2/sched/ 
# Drudg vex file
ssh oper@pcfs-2hb 'cd /usr2/sched; echo -e "hb\n11\n19 16 1 1\n\n12\n0\n"|/usr2/fs/bin/drudg "$1.vex /usr2/fs/bin/drudg'
# Copy template proc file to /usr2/proc
ssh oper@pcfs-2hb 'cp /usr2/proc/spirals.prc /usr2/proc/$1hb.prc'

echo "Katherine"
# Copy vex to pcfske home area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfske:~/
# drudg the vex file there
ssh oper@pcfshb 'cd /usr2/sched; echo -e "ke\n3\n0\n"|/usr2/fs/bin/drudg "$1.vex'
# Run the snp2flex script
ssh oper@pcfshb snp2flex.py $1ke.snp

# Copy vex file over the pcfs-2ke sched area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2ke:/usr2/sched/ 
# Drudg vex file
ssh oper@pcfs-2ke 'cd /usr2/sched; echo -e "ke\n11\n19 16 1 1\n\n12\n0\n"|/usr2/fs/bin/drudg "$1.vex /usr2/fs/bin/drudg'
# Copy template proc file to /usr2/proc
ssh oper@pcfs-2ke 'cp /usr2/proc/spirals.prc /usr2/proc/$1ke.prc'

echo "Setting up correlator"
cp -r /home/correlations/template /home/correlations2/$1/
cp /home/observer/correlations2/spiralssched/vex/$1.vex /home/observer/correlations2/$1/$1.vex
cd /home/observer/correlations2/$1/ && ../../corr_scripts/kludge_vex.sh $1

touch hbX.filelist
touch hbY.filelist
touch keX.filelist
touch keY.filelist
touch cd.filelist
touch wa.filelist
touch ho.filelist
touch wa.filelist
touch yg.filelist

echo "Done"
