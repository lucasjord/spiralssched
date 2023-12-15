#!/bin/bash

cd /home/observer/correlations2/spiralssched/ && git pull

# Check date and time of experiment match current times

export START=$(grep exper_nominal_start /home/observer/correlations2/spiralssched/year3/vex/$1.vex | cut -d '=' -f 2 | tr -d ";")
export STOP=$(grep exper_nominal_stop /home/observer/correlations2/spiralssched/year3/vex/$1.vex | cut -d '=' -f 2 | tr -d ";")

export YEAR=$(echo $START | cut -d 'y' -f 1)
export DOY=$(echo $START | cut -d 'y' -f 2 | cut -d 'd' -f 1)

export ffsource=$(grep source= /home/observer/correlations2/spiralssched/year3/vex/$1.vex |cut -d ';' -f3 |head -1 |cut -d '=' -f2)

# Check year, exit if no match
if [ $YEAR -eq $(date -u +%Y) ]; then :
else echo 'Experiment does not take place in '$(date -u +%Y%j%H%M%S); exit
fi

# Check DOY, requires response if incorrect
if [ $DOY -eq $(date -u +%j) ]
	then :
else 
    echo 'Experiment '$1' does not take place today:'
    echo "Experiment start: $START"
    echo "Time is now:      $(date -u +%Yy%jd%Hh%Mm%Ss)"
    echo "Is this okay? [y|n]'"
	read RESP
	while [ ${#RESP[@]} -eq 0 ];
		do echo 'Please answer with y or n'
		read RESP
	done
   	if [[ "$RESP" == y* ]] || [[ "$RESP" == Y* ]]
      		then :
      	else echo 'Exiting, please check vex date.'; exit
   	fi
fi

# Identify experiment series/class
export class=$(echo $1 | sed 's/.$//')

array=$(grep Array /home/observer/correlations2/spiralssched/year3/vex/$1.vex | cut -d ":" -f 3)
echo $array


if [[ *$array* == *"Cd"* ]]; then
    echo ' '
    echo "############################################################"
    echo "Ceduna"
    echo ' '
    # Copying vex file to Ceduna
    scp /home/observer/correlations2/spiralssched/year3/vex/$1.vex oper@pcfscd:/usr2/sched/$1.vex
    # Removing old files
    ssh oper@pcfscd 'cd /usr2/sched; rm '$1cd.*
    # Drudging
    ssh oper@pcfscd 'cd /usr2/sched; echo -e "cd\n11\n19 16 1 1\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
    # Copying template vex over
    ssh oper@pcfscd 'cp /usr2/proc/spirals.prc /usr2/proc/'$1cd.prc
fi

if [[ *$array* == *"Hb"* ]]; then
    echo ' '
    echo "############################################################"
    echo "Hobart12"
    echo ' '
    # Copy vex to pcfshb sched area
    scp /home/observer/correlations2/spiralssched/year3/vex/$1.vex oper@pcfshb:/usr2/sched/
    # drudg the vex file there
    ssh oper@pcfshb rm /usr2/sched/$1hb.snp
    ssh oper@pcfshb 'cd /usr2/sched; echo -e "hb\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
    # Copy template proc file over
    ssh oper@pcfshb 'cp /usr2/proc/spirals.prc /usr2/proc/'$1hb.prc
fi

if [[ *$array* == *"Ke"* ]]; then
    echo ' '
    echo "############################################################"
    echo "Katherine"
    echo ' '
    # Copy vex to pcfske sched area
    scp /home/observer/correlations2/spiralssched/year3/vex/$1.vex oper@pcfske:/usr2/sched/
    # drudg the vex file there
    ssh oper@pcfske 'cd /usr2/sched; rm '$1ke.snp
    ssh oper@pcfske 'cd /usr2/sched; echo -e "ke\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
    ssh oper@pcfske 'cp /usr2/proc/spirals.prc /usr2/proc/'$1ke.prc
fi

if [[ *$array* == *"Yg"* ]]; then
    echo ' '
    echo "############################################################"
    echo "Yarragadee"
    echo ' '
    # Copy vex to pcfske sched area
    scp /home/observer/correlations2/spiralssched/year3/vex/$1.vex oper@pcfs-2yg:/usr2/sched/
    # drudg the vex file there
    ssh oper@pcfs-2yg 'cd /usr2/sched; rm '$1yg.snp
    ssh oper@pcfs-2yg 'cd /usr2/sched; echo -e "yg\n3\n0\n"|/usr2/fs/bin/drudg '$1.vex
    ssh oper@pcfs-2yg 'cp /usr2/proc/spirals.prc /usr2/proc/'$1yg.prc
fi

echo ' '
echo "############################################################"
echo "Setting up correlator"
echo ' '
mkdir -p /home/observer/correlations2/$1
#cp  /home/observer/correlations2/template/threads /home/observer/correlations2/$1/threads
#cp  /home/observer/correlations2/template/machines /home/observer/correlations2/$1/machines
cp  /home/observer/correlations2/template/template2.v2d /home/observer/correlations2/$1/$1.v2d
cp /home/observer/correlations2/spiralssched/year3/vex/$1.vex /home/observer/correlations2/$1/$1.vex
cd /home/observer/correlations2/$1/ && /home/observer/correlations2/corr_scripts/kludge_vex.sh $1

echo "Making v2d file"
/home/observer/correlations2/corr_scripts/eop.sh $YEAR $DOY >> /home/observer/correlations2/$1/$1.v2d
sed -i "s/FAKESTART/$START/g" /home/observer/correlations2/$1/$1.v2d
sed -i "s/FAKESTOP/$STOP/g" /home/observer/correlations2/$1/$1.v2d
sed -i "s/FAKEVEX/$1/g" /home/observer/correlations2/$1/$1.v2d
sed -i "s/FAKESOURCE/$ffsource/g" /home/observer/correlations2/$1/$1.v2d

# Making empty filelists
touch /home/observer/correlations2/$1/hbX.filelist
touch /home/observer/correlations2/$1/hbY.filelist
touch /home/observer/correlations2/$1/keX.filelist
touch /home/observer/correlations2/$1/keY.filelist
touch /home/observer/correlations2/$1/cd.filelist
touch /home/observer/correlations2/$1/wa.filelist
touch /home/observer/correlations2/$1/ho.filelist
touch /home/observer/correlations2/$1/ygX.filelist
touch /home/observer/correlations2/$1/ygY.filelist

echo ' '
echo "##############################"
echo "Done"
echo ' '
echo "Experiment starts $START"
echo "Time is now       $(date -u +%Yy%jd%Hh%Mm%Ss)"
echo ' '

