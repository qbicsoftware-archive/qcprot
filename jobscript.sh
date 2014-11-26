#!/bin/bash
# properties = {properties}

set -e

WORKDIR=$(readlink -f $1)
# Directory where the script is stored.
# See https://stackoverflow.com/questions/59895/
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


if [ ! -f ${SCRIPTDIR}/Snakefile ] ; then
	echo "Could not find snakefile in ${SCRIPTDIR}" 1>&2
	exit 1
fi


module load lib/openms/1.11
module load bio/xtandem
module load bio/myrimatch
module load math/R

# Anaconda is a python distribution (http://continuum.io/)
# Snakemake must be installed.
export PATH=/lustre_cfc/software/qbic/anaconda3/bin:$PATH

export R_HOME=${SCRIPTDIR}/r_scripts
export QCPROT_VERSION=$(git -C ${SCRIPTDIR} describe --always --dirty)
export OPENMS_BIN=$(dirname "$(which IDMapper)")

{exec_job}
