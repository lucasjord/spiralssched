#!/bin/bash

cd /home/observer/correlations2/spiralssched/ && git pull

# Check date and time of experiment match current times

export START=$(grep exper_nominal_start /home/observer/correlations2/spiralssched/vex/$1.vex | cut -d '=' -f 2 | tr -d ";")
export STOP=$(grep exper_nominal_stop /home/observer/correlations2/spiralssched/vex/$1.vex | cut -d '=' -f 2 | tr -d ";")

export YEAR=$(echo $START | cut -d 'y' -f 1)
export DOY=$(echo $START | cut -d 'y' -f 2 | cut -d 'd' -f 1)

# Check year, exit if no match
if [ $YEAR -eq $(date -uj +%Y) ]; 
	:
else 
	echo 'Experiment does not take place in '$(date -uj +%Y); exit
fi

# Check DOY, requires response if incorrect
if [ $DOY -eq $(date -uj +%j) ]
	:
else 
	then 
		echo 'Experiment '$1' does not take place today: '$DOY' vs. '$(date -uj +%j)', is this okay? [y|n]'
		read RESP
		while [ ${#RESP[@]} -eq 0 ]; 
			do echo 'Please answer with y or n'
			read RESP
		done
   		if [[ "$RESP" == y* ]] || [[ "$RESP" == Y* ]]
      		:
      	else 
      		echo 'Exiting, please check vex date.'; exit
   		fi 
fi

# Identify experiment series/class
export class=$(echo $1 | sed 's/.$//')

echo ' '
echo "############################################################"
echo "Ceduna"
echo ' '
# Copying vex file to Ceduna
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfscd:/usr2/sched/$1.vex
# Removing old files
ssh oper@pcfscd 'cd /usr2/sched; rm '$1cd.*
# Drudging
ssh oper@pcfscd 'cd /usr2/sched; echo -e "cd\n11\n19 16 1 1\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
# Add source=stow
ssh oper@pcfscd 'echo "source=stow" >> /usr2/sched/'$1cd.snp
# Copying template vex over
ssh oper@pcfscd 'cp /usr2/proc/spiralscd.prc /usr2/proc/'$1cd.prc

echo ' '
echo "############################################################"
echo "Hobart12"
echo ' '
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
ssh oper@pcfs-2hb 'cd /usr2/sched; echo -e "hb\n11\n19 16 1 1\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
# Add source=stow
ssh oper@pcfshb 'echo "source=stow" >> /usr2/sched/'$1hb.snp
# Copy template proc file to /usr2/proc
ssh oper@pcfs-2hb 'cp /usr2/proc/spirals.prc /usr2/proc/'$1hb.prc

echo ' '
echo "############################################################"
echo "Katherine"
echo ' '
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
ssh oper@pcfs-2ke 'cd /usr2/sched; echo -e "ke\n11\n19 16 1 1\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
# Add source=stow
ssh oper@pcfs-2ke 'echo "source=stow" >> /usr2/sched/'$1ke.snp
# Copy template proc file to /usr2/proc
ssh oper@pcfs-2ke 'cp /usr2/proc/spirals.prc /usr2/proc/'$1ke.prc

echo ' '
echo "############################################################"
echo "Setting up correlator"
echo ' '
mkdir -p /home/observer/correlations2/$1
cp /home/observer/correlations2/template/* /home/observer/correlations2/$1/*
cp /home/observer/correlations2/spiralssched/vex/$1.vex /home/observer/correlations2/$1/$1.vex
/home/observer/correlations2/corr_scripts/kludge_vex.sh $1

# Making v2d file
mv /home/observer/correlations2/$1/template.v2d /home/observer/correlations2/$1/$1.v2d
/home/observer/correlations2/corr_scripts/eop.sh $YEAR $DOY >> /home/observer/correlations2/$1/$1.v2d
sed -i "s/FAKESTART/$START/g" /home/observer/correlations2/$1/$1.v2d
sed -i "s/FAKESTOP/$STOP/g" /home/observer/correlations2/$1/$1.v2d
sed -i "s/FAKEVEX/$1/g" /home/observer/correlations2/$1/$1.v2d

# Making empty filelists
touch /home/observer/correlations2/$1/hbX.filelist
touch /home/observer/correlations2/$1/hbY.filelist
touch /home/observer/correlations2/$1/keX.filelist
touch /home/observer/correlations2/$1/keY.filelist
touch /home/observer/correlations2/$1/cd.filelist
touch /home/observer/correlations2/$1/wa.filelist
touch /home/observer/correlations2/$1/ho.filelist
touch /home/observer/correlations2/$1/wa.filelist
touch /home/observer/correlations2/$1/yg.filelist

echo ' '
echo "##############################"
echo "Done"
echo ' '
