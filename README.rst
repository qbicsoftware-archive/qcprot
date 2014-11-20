Run a quality control workflow on mzML files from mass spec based proteomics.

Execute `qcprot` with the working directory as parameter. This directory must
contain

- mzML : A directory that contains the mzML files
- config.json : a json file with the parameters for the qc workflow.

Results will be written to a directory `results` in the workdir.

Excample config file:

    {
        "REF" : "uniprot_sprot_101104_human_concat.fasta",
        "input_files" : ["mzML/velos005614.mzML"],
        "precursor_tolerance" : "15",
        "precursor_tolerance_unit" : "ppm",
        "fragment_tolerance" : "25",
        "fragment_tolerance_unit" : "ppm",
        "retention_time_tolerance": "60",
        "retention_time_tolerance_unit" : "sec",
        "shell_prefix" : "module load lib/openms/1.11"
    }
