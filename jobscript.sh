#!/bin/sh
#PBS -q cfc
#PBS -A qbic
#PBS -l nodes=1:ppn=10:cfc
#PBS -l walltime=40:00:00
# properties = {properties}

source /lustre_cfc/ws_qstor/ws/qeaco01-conda_openms-0/usr/bin/activate /lustre_cfc/ws_qstor/ws/qeaco01-conda_openms-0/usr
module load lib/openms/1.11
module load bio/xtandem/201309011-gnu-4.8
module load math/R/3.1.0

{exec_job}
exit 0
