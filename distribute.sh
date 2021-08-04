#!/bin/bash

cd /home/observer/correlations2/spiralssched/ && git pull

export class=$(echo $1 | sed 's/.$//')

echo ' '
echo "##############################"
echo "Newsmerd for slogit"
# Copying vex file to newsmerd for slogit
ssh observer@newsmerd mkdir /vlbobs/lba/$class
scp /home/observer/correlations2/spiralssched/vex/$1.vex observer@newsmerd.phys.utas.edu.au:/vlbobs/lba/$class/$1.vex

echo "##############################"
echo "Hobart12"
# Copy vex to pcfshb home area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfshb:~/
# drudg the vex file there
ssh oper@pcfshb rm $1hb.snp
ssh oper@pcfshb 'echo -e "hb\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
# Run the snp2flex script
ssh oper@pcfshb /usr2/oper/AuscopeVLBI/fs/src/snp2flexbuff.py $1hb.snp

# Copy vex file over the pcfs-2hb sched area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2hb:/usr2/sched/ 
# Drudg vex file
ssh oper@pcfs-2hb 'cd /usr2/sched/; rm '$1hb.snp
ssh oper@pcfs-2hb 'cd /usr2/sched; echo -e "hb\n11\n19 16 1 1\n3\n\n0\n"|/usr2/fs/bin/drudg '$1.vex
# Copy template proc file to /usr2/proc
ssh oper@pcfs-2hb 'cp /usr2/proc/spirals.prc /usr2/proc/'$1hb.prc

echo ' '
echo "##############################"
echo "Katherine"
# Copy vex to pcfske home area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfske:~/
# drudg the vex file there
ssh oper@pcfske rm $1ke.snp
ssh oper@pcfske 'echo -e "ke\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
# Run the snp2flex script
ssh oper@pcfske /usr2/oper/AuscopeVLBI/fs/src/snp2flexbuff.py $1ke.snp

# Copy vex file over the pcfs-2ke sched area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2ke:/usr2/sched/ 
# Drudg vex file
ssh oper@pcfs-2ke 'cd /usr2/sched; rm '$1ke.snp
ssh oper@pcfs-2ke 'cd /usr2/sched; echo -e "ke\n11\n19 16 1 1\n\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
# Copy template proc file to /usr2/proc
ssh oper@pcfs-2ke 'cp /usr2/proc/spirals.prc /usr2/proc/'$1ke.prc

echo ' '
echo "##############################"
echo "Setting up correlator"
cp -r /home/observer/correlations2/template /home/observer/correlations2/$1/
cp /home/observer/correlations2/spiralssched/vex/$1.vex /home/observer/correlations2/$1/$1.vex
cd /home/observer/correlations2/$1/ && /home/observer/correlations2/corr_scripts/kludge_vex.sh $1

# Making empty filelists
touch hbX.filelist
touch hbY.filelist
touch keX.filelist
touch keY.filelist
touch cd.filelist
touch wa.filelist
touch ho.filelist
touch wa.filelist
touch yg.filelist

echo ' '
echo "##############################"
echo "Done"
