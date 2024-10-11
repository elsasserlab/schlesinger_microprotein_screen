#!/bin/bash -l

#SBATCH -A naiss2023-22-84
#SBATCH -p devel
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -M rackham

module load mono

mono /proj/snic2020-6-3/MaxQuant/bin/MaxQuantCmd.exe /proj/snic2020-6-3/Doerte/MS_DS_231027/20231027_6reps.xml
