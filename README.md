# motif

Identify motifs for CLIP peak clusters.

# Usage
1. Export this path `/storage/vannostrand/software/motif/venv/bin` to your `PATH`.
```shell script
$ export PATH=/storage/vannostrand/software/motif/venv/bin:$PATH
```
2. Check the usage:
```shell script
$ motif -h
usage: motif [-h] BED SPECIES OUTDIR UID

Identify motifs for CLIP peak clusters.

positional arguments:
  BED         Path to the BED file, usually this is the peak clusters normalized compressed bed file.
  SPECIES     Species shortname, can be one of these: hg19, hg39, mm10.
  OUTDIR      Path to the output directory, default: the BED file's parent directory.
  UID         User specified identifier for the dataset.

optional arguments:
  -h, --help  show this help message and exit
```

3. Pass required arguments to `motif` to start motif identification.

