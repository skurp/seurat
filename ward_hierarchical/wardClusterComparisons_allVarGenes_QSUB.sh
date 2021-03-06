#!/bin/bash

# Will Connell
# 2017-12-15
# Code to run wardClusterComparisons_allVarGenes.R

# Reminder:
# cd dkohn/wconnell/geschwind/projects/SC3/Scripts
#  make /logs directory in code directory
# Must load modules and qsub from login node:
#  module load gcc/4.9.3
#  module load R/3.4.0

# qsub:
# qsub wardClusterComparisons_allVarGenes_QSUB.sh
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -N wardClusterComparisons_allVarGenes
#$ -o logs/wardClusterComparisons_allVarGenes_QSUB_$JOB_ID_$TASK_ID.log
#$ -e logs/wardClusterComparisons_allVarGenes_QSUB_$JOB_ID_$TASK_ID.error
#$ -l h_data=8,h_rt=4:00:00 -pe shared 8
#$ -t 1-1
#$ -m bea
################################################################################
echo ""
echo "Starting wardClusterComparisons_allVarGenes_QSUB.sh ${SGE_TASK_ID}... "$(date)
echo ""
##########################################################################

# Path to R
# pathRscript=/u/local/apps/R/3.3.0/gcc-4.4.7/enable-R-shlib/bin/Rscript
pathRscript=mod/R/3.4.0/gcc-4.9.3_MKL-2017.0/lib64/R/bin/Rscript
# Set path to GCC for Seurat package
PATH=/u/local/compilers/gcc/4.9.3/bin:$PATH
INCLUDE=/u/local/compilers/gcc/4.9.3/include:$INCLUDE
LD_LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LD_LIBRARY_PATH
LIBRARY_PATH=/u/local/compilers/gcc/4.9.3/lib:/u/local/compilers/gcc/4.9.3/lib64:$LIBRARY_PATH

## Run TF_Expression.R
${pathRscript} wardClusterComparisons_allVarGenes.R
##########################################################################

echo ""
echo "40k_mergedInfomapLJ_25DEGenes_QSUB.sh ${SGE_TASK_ID}... "$(date)
##########################################################################
