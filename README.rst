Run a quality control workflow on mzML files from mass spec based proteomics.

Execute ``qcprot`` with the working directory as parameter. This directory must
contain a config file called ``config.json`` and a directory ``mzml`` that
contains the input mzml files.

Results will be written to a directory ``results`` in the workdir.

Example config file::

    {
        "fasta" : "uniprot_sprot_101104_human_concat.fasta",
        "fasta_has_decoy": false,
        "ini_path" : "path/to/ini/files"
    }

Only ``"fasta"`` is mandatory. If you want to include several fasta files
(``'crap.fasta'!``) just specify a list of fasta files::

    {
        "fasta" : ["uniprot_sprot_101104_human_concat.fasta", "crap.fasta"]
    }

QBiC specific
-------------

Load ``qcprot`` on the cluster with::

    module load qbic/qcprot

Create a working directory as described above and execute the following::

    qcprot /path/to/workdir [Snakemake args]

For example::

    # Execute up to 2 task at the same time (*not* on the frontend!)
    qcprot /path/to/workdir -j2

    # Use up to 10 nodes of the cluster (execute on the frontend)
    qcprot workdir --cluster "qsub -A qbic -l nodes=1:ppn=20:cfc -q cfc -V" -j10
