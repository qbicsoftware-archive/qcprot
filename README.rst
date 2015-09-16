Run a quality control workflow on mzML files from mass spec based proteomics.

This workflow is best executed through `qproject`.
(https://github.com/qbicsoftware/qproject)

Create the directory structure of the workdir (it should not exist):

```
qproject create -t path/to/workdir -w github:qbicsoftware/qcprot
```

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
