#!/bin/bash

# Will Connell
# 2017-12-14
# Code to run seuratDatasetAlignment_full_dataset.R

# Reminder:
# cd dkohn/wconnell/geschwind/projects/SC3/Scripts
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.4.0

# qsub:
# qsub seuratDatasetAlignment_full_dataset_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N seuratDatasetAlignment_full_dataset
#$ -o logs/seuratDatasetAlignment_full_dataset_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/seuratDatasetAlignment_full_dataset_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=8G,h_rt=8:00:00 -pe shared 6
#$ -t 1-1
#$ -m bea
################################################################################
echo ""
echo "Starting seuratDatasetAlignment_full_dataset_QSUB.sh ${SGE_TASK_ID}... "$(date)
echo ""
##########################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=/u/local/apps/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

## Run seuratDatasetAlignment_full_dataset.R
${pathRscript} seuratDatasetAlignment_full_dataset.R
##########################################################################

echo ""
echo "seuratDatasetAlignment_full_dataset_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
