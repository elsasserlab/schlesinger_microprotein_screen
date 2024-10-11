#!/bin/bash -l

#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 8-1:00:00
#SBATCH -M Rackham

module load mono

mono /proj/snic2020-6-3/MaxQuant/bin/MaxQuantCmd.exe /proj/snic2020-6-3/Doerte/MS_DS_231027/20231027_6reps.xml
