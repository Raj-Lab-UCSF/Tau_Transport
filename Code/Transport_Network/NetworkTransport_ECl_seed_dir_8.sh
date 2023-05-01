#!/bin/bash
#### Specify job name
#SBATCH -J Tau_Transport
#### Output file
#SBATCH -o "%x"_"%j".out
#### Error file
#SBATCH -e "%x"_"%j".err
#### number of cores 
#SBATCH -n 16
#### Specify queue
#SBATCH --partition=long
#### --nodelist=oakland,piedmont,novato,quartzhill
#### memory per core
#SBATCH --mem=8G
#### Maximum run time 
#SBATCH --time=2-00:00:00

matlab -nosplash -nodesktop -r "run('~/Tau_Transport/Code/Transport_Network/NetworkTransportModel_Wrapper_GridSearch_ECl_seed_8.m'); quit"
