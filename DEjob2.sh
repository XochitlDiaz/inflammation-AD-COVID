#!/bin/bash
#
#
# Job Name
#PBS -N DE.job
#
# Resources
#PBS -l ncpus=2
#PBS -l mem=36gb
#PBS -l walltime=8:00:00
# Send an email after the jon has finished
#PBS -m e
#PBS -M xochitl.happy@gmail.com
#
# Load modules
module purge
module load R/4.2.0
#
cd $PBS_O_WORKDIR
# Write commands 
R -f DE_acrossctrlAD.R
