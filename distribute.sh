#!/bin/bash

cd /home/observer/correlations2/spiralssched/ && git pull 

export class=$(echo $1 | sed 's/.$//')

echo "\n Newsmerd for slogit"
rsync -Puv /home/observer/correlations2/spiralssched/vex/$1.vex observer@newsmerd:/vlbobs/lba/$class/$1.vex
echo "\n Hobart12"
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfshb:~/
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2hb:/usr2/sched/
echo "\n Katherine"
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfske:~/
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2ke:/usr2/sched/

echo "\n hex6"
rsync -Puv /home/observer/correlations2/spiralssched/vex/$1.vex /home/observer/correlations2/$1/$1.vex

echo "\n Done"
