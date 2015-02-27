#!/bin/sh
# properties = {properties}

module load devel/qt/4.8.4
module load qbic/openms/1.11.1-2996-g13ffbd7
module load bio/xtandem/201309011-gnu-4.8
module load math/R/3.1.0

{exec_job}
exit 0
