#! /bin/bash

#/Applications/MATLAB_R2024b.app/bin/matlab -r "forward_model_cluster 0.000001 0.0005 0.05 20 20 1"

/Applications/MATLAB_R2024b.app/bin/matlab -r "forward_model_cluster $1 $2 $3 $4 $5 $6"

#module load matlab
#matlab -r "forward_model_cluster $1 $2 $3 $4 $5 $6"