Run a quality control workflow on mzML files from mass spec based proteomics.

Execute ``qcprot`` with the working directory as parameter. This directory must
contain a config file called ``config.json`` and a directory ``mzml`` that
contains the input mzml files.

Results will be written to a directory ``results`` in the workdir.

Example config file::

    {
        "fasta" : "uniprot_sprot_101104_human_concat.fasta",
        "fasta_has_decoy": false,
        "ini_path" : "path/to/ini/files",
    }

Only ``"fasta"`` is mandatory.
