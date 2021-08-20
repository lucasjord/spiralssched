#!/bin/bash

cd /home/observer/correlations2/spiralssched/ && git pull

# Check date and time of experiment match current times

export START=$(grep exper_nominal_start /home/observer/correlations2/spiralssched/vex/$1.vex | cut -d '=' -f 2 | tr -d ";")
export STOP=$(grep exper_nominal_stop /home/observer/correlations2/spiralssched/vex/$1.vex | cut -d '=' -f 2 | tr -d ";")

export YEAR=$(echo $START | cut -d 'y' -f 1)
export DOY=$(echo $START | cut -d 'y' -f 2 | cut -d 'd' -f 1)

# Check year, exit if no match
if [ $YEAR -eq $(date -u +%Y) ]; then :
else echo 'Experiment does not take place in '$(date -uj +%Y); exit
fi

# Check DOY, requires response if incorrect
if [ $DOY -eq $(date -u +%j) ]
	then :
else echo 'Experiment '$1' does not take place today: '$DOY' vs. '$(date -u +%j)', is this okay? [y|n]'
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

echo ' '
echo "############################################################"
echo "Ceduna"
echo ' '
# Copying vex file to Ceduna
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfscd:/usr2/sched/$1.vex
ssh oper@pcfscd << EOF
	# Removing old files
	cd /usr2/sched; rm $1cd.*
	# Drudging
	cd /usr2/sched; echo -e "cd\n11\n19 16 1 1\n3\n0\n"|/usr2/fs/bin/drudg $1.vex
	# Add source=stow
	echo "source=stow" >> /usr2/sched/$1cd.snp
	# Copying template vex over
	cp /usr2/proc/spiralscd.prc /usr2/proc/$1cd.prc
EOF

echo ' '
echo "############################################################"
echo "Hobart12"
echo ' '
# Copy vex to pcfshb home area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfshb:~/
ssh oper@pcfhb << EOF
	# drudg the vex file there
	rm $1hb.snp
	echo -e "hb\n3\n0\n"|/usr2/fs/bin/drudg $1.vex
	# Run the snp2flex script
	/usr2/oper/AuscopeVLBI/fs/src/snp2flexbuff.py $1hb.snp
EOF

# Copy vex file over the pcfs-2hb sched area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2hb:/usr2/sched/ 
# Drudg vex file
ssh oper@pcf-2hb << EOF
	rm /usr2/sched/$1hb.snp
	cd /usr2/sched; echo -e "hb\n11\n19 16 1 1\n3\n0\n"|/usr2/fs/bin/drudg $1.vex
	# Add source=stow
	echo "source=stow" >> /usr2/sched/$1hb.snp
	# Copy template proc file to /usr2/proc
	cp /usr2/proc/spirals.prc /usr2/proc/$1hb.prc
EOF

echo ' '
echo "############################################################"
echo "Katherine"
echo ' '

# Copy vex to pcfske home area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfske:~/
ssh oper@pcfke << EOF
	# drudg the vex file there
	rm $1ke.snp
	echo -e "ke\n3\n0\n"|/usr2/fs/bin/drudg $1.vex
	# Run the snp2flex script
	/usr2/oper/AuscopeVLBI/fs/src/snp2flexbuff.py $1ke.snp
EOF

# Copy vex file over the pcfs-2ke sched area
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2ke:/usr2/sched/ 
# Drudg vex file
ssh oper@pcfs-2ke << EOF
	rm /usr2/sched/$1ke.snp
	cd /usr2/sched; echo -e "ke\n11\n19 16 1 1\n3\n0\n"|/usr2/fs/bin/drudg $1.vex
	# Add source=stow
	echo "source=stow" >> /usr2/sched/$1ke.snp
	# Copy template proc file to /usr2/proc
	cp /usr2/proc/spirals.prc /usr2/proc/$1ke.prc
EOF

echo ' '
echo "############################################################"
echo "Setting up correlator"
echo ' '
mkdir -p /home/observer/correlations2/$1
cp /home/observer/correlations2/template/threads /home/observer/correlations2/$1/threads
cp /home/observer/correlations2/template/machines /home/observer/correlations2/$1/machines
cp /home/observer/correlations2/template/template.v2d /home/observer/correlations2/$1/$1.v2d
cp /home/observer/correlations2/spiralssched/vex/$1.vex /home/observer/correlations2/$1/$1.vex
cd /home/observer/correlations2/$1/ && /home/observer/correlations2/corr_scripts/kludge_vex.sh $1

echo "Making v2d file"
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
echo "Experiment starts $START"
echo "Time is now $(date -u +%Y.%j.%H.%M.%S)"
echo ' '
