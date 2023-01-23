#!/bin/bash
#PBS -q himem
#PBS -l nodes=1:ppn=16
#PBS -l mem=240gb
#PBS -l walltime=300:00:00

#!PBS -l nodes=1:ppn=16,walltime=100:00:00

#! change the working directory (default is home directory)
cd $PBS_O_WORKDIR
#! Record some useful job details in the output file
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo PBS job ID is $PBS_JOBID
echo This jobs runs on the following nodes:
echo `cat $PBS_NODEFILE | uniq`
#! add the MATLAB module (as per BCp3)
module add apps/matlab-r2016b
options="-nodesktop -nosplash -noFigureWindows"
#! Run MATLAB in batch mode
matlab $options -r BatchMicAnalysisTaskOnly_Unsorted	