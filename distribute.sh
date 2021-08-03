#!/bin/bash

git pull 

export class=$(sed 's/.$//' $1)
rsync -Puv /home/observer/correlations2/spiralssched/vex/$1.vex observer@newsmerd:/vlbobs/lba/$class/$1.vex
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfshb:~/
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2hb:/usr2/sched/
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfske:~/
scp /home/observer/correlations2/spiralssched/vex/$1.vex oper@pcfs-2ke:~/

