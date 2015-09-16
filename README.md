Run a quality control workflow on mzML files from mass spec based proteomics.

## Requirements
- A python 3 environment with (eg. anaconda)
  - snakemake
  - jinja2
  - qproject (eg. `pip install git+https://github.com/qbicsoftware/qproject`)
  - xtandem
  - R
- openms 2.0 (there is a conda package that can be installed with `conda
  install -c aseyboldt openms` that might or might not work)

##Usage
Create the directory structure of the workdir (it should not exist):

```
qproject create -t path/to/workdir -w github:qbicsoftware/qcprot
```

This will create the `<workdir>` and some directories inside and clone
`qcprot` to `<workdir>/src`.
Copy the mzML files to `<workdir>/data` and your fasta file(s) to
`<workdir>/ref`.

Copy the ini files to `etc`:

```
cd <workdir>/src/inis <workdir>/etc
```

and modify if necessary.

Execute the workflow with

```
cd <workdir>/src
snakemake
```

You can also adjust the jobscript and use qproject to execute the workflow:

```
qproject run -t <workdir>
```

If qproject is used to execute `qcprot` some data about the run is
stored in `<workdir>/archive` for reproducibility. This is still
very much work in progress though' The output of `snakemake`
will end up in `<workdir>/logs/snake.err`.
