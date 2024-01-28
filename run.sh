#!/bin/bash

#SBATCH -p cpu3
#SBATCH -N 1-1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -o 1.out
#SBATCH --nodelist wmc-slave-g12

date
module load matlab


srun matlab -nodisplay -nosplash -nodesktop -r "clear; \
stego_dir = 'HILL_payload0.5_lambda0.15_N25_coef13'; \
lambda=0.15; \
coef = 13; \
payload = 0.5; \
mean_dir = '/data-x/g17/zhangjiansong/25/mean'; \
var_dir = '/data-x/g17/zhangjiansong/25/var'; \
run('example.m'); exit;"









