# motif

Identify motifs for CLIP peak clusters.

# Usage
1. On TACO cluster, directory call `/storage/vannostrand/software/motif/motif` 
with one required parameter `annotated_peak_bed`, a path to a 10-column bed file 
contains annotated peaks to identify motifs for peak clusters. In case you need 
to fine tune the results, pass any of the 4 optional parameters (See check the 
usage section below for details).


2. Check the usage:
```shell script
$ /storage/vannostrand/software/motif/motif

Identify motifs for eCLIP peak clusters.

Usage:
motif annotated_peak_bed [species] [outdir] [uid] [l10p] [l2fc] [processes]
    - annotated_peak_bed: 10-column bed file contains annotated peaks.
    - species: species name, e.g., hg19, hg38, ..., default: hg19.
    - outdir: path to output directory, default: current work directory.
    - uid: unique identifier for the dataset,
           default: annotated_peak_bed file basename without .bed extension.
    - l10p: l10p cutoff, default: 3.
    - l2fc: l2fc cutoff, default: 3.
    - processes: Number of processes to use, default: 1.
```

## Example
1. Check annotated peak bed file. Assume we have an annotated peak bed file: 
`/storage/vannostrand/analysis/20210908_Eric/DHX16_seCLIP/DHX16.F1.ip.peak.clusters.normalized.compressed.annotated.bed`, 
which is a 10-column bed file contains annotated peaks:

```
$ cd /storage/vannostrand/analysis/20210908_Eric/DHX16_seCLIP
$ head -n 5 DHX16.F1.ip.peak.clusters.normalized.compressed.annotated.bed
chrM    2019    2027    400     6.08911863559215        +       noncoding_exon;ENSG00000210082.2        noncoding_exon||ENSG00000210082.2       ENSG00000210082.2       MT-RNR2 noncoding_exon|8
chrM    2013    2019    400     5.91884366567356        +       noncoding_exon;ENSG00000210082.2        noncoding_exon||ENSG00000210082.2       ENSG00000210082.2       MT-RNR2 noncoding_exon|6
chrM    2027    2031    400     5.89435586261036        +       noncoding_exon;ENSG00000210082.2        noncoding_exon||ENSG00000210082.2       ENSG00000210082.2       MT-RNR2 noncoding_exon|4
chrM    1971    2000    400     5.89069384062431        +       noncoding_exon;ENSG00000210082.2        noncoding_exon||ENSG00000210082.2       ENSG00000210082.2       MT-RNR2 noncoding_exon|29
chrM    2000    2013    400     5.82421124135645        +       noncoding_exon;ENSG00000210082.2        noncoding_exon||ENSG00000210082.2       ENSG00000210082.2       MT-RNR2 noncoding_exon|13
```

2. Call `motif` to identify motifs.
Inside `/storage/vannostrand/analysis/20210908_Eric/DHX16_seCLIP`, we can call `motif` as follows to identify 
motifs and save results to a current directory (`$PWD`) and dataset set to `DHX16.F1`:
```shell
$ /storage/vannostrand/software/motif/motif \
  DHX16.F1.ip.peak.clusters.normalized.compressed.annotated.bed \
  hg19 \
  $PWD \
  DHX16.F1
```

3. After successfully run the above command, you should be able to see a directory named 
`DHX16.F1.motifs.40.min.plus.5p.20.3p.5` inside the current directory and it has all the 
results for identified motifs:
```shell
$ ls DHX16.F1.motifs.40.min.plus.5p.20.3p.5
3utr
5ss
5utr
all
CDS
distintron
miRNA
miRNA_proximal
noncoding_distintron
noncoding_exon
proxintron
tRNA
DHX16.F1.3utr.homer.results.html
DHX16.F1.5ss.homer.results.html
DHX16.F1.5utr.homer.results.html
DHX16.F1.all.homer.results.html
DHX16.F1.CDS.homer.results.html
DHX16.F1.distintron.homer.results.html
DHX16.F1.miRNA.homer.results.html
DHX16.F1.miRNA_proximal.homer.results.html
DHX16.F1.noncoding_distintron.homer.results.html
DHX16.F1.noncoding_exon.homer.results.html
DHX16.F1.proxintron.homer.results.html
DHX16.F1.tRNA.homer.results.html
```


