#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#### Specify job name
#$ -N Tau_Transport
#### Output file
#$ -o /wynton/protected/home/rajlab/jtorok/JobOutputs/$JOB_NAME_$JOB_ID.out
#### Error file
#$ -e /wynton/protected/home/rajlab/jtorok/JobOutputs/$JOB_NAME_$JOB_ID.out
#### number of cores
#$ -pe smp 32
#### memory per core
#$ -l mem_free=2G
#### Maximum run time
#$ -l h_rt=150:00:00

module load matlab
matlab -batch "NetworkTransport_ECl_hippocampome_final_round2_2"