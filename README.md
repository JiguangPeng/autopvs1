# AutoPVS1
An automatic classification tool for PVS1 interpretation of null variants.

A web version for AutoPVS1 is also provided: http://autopvs1.genetics.bgi.com/

## PREREQUISITE
### 1. Variant Effect Predictor (VEP)
**autopvs1** use [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) to determine the effect of your 
variants (SNVs, insertions, deletions, CNVs) on genes, 
transcripts, and protein sequence.

To get HGVS name for the variant, homo_sapiens_refseq 97_GRCh37 cache should be used.

```bash
# VEP Installation
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git pull
git checkout release/97
perl INSTALL.pl
```

### 2. pyfaidx
Samtools provides a function “faidx” (FAsta InDeX), which creates a small flat index file “.fai” 
allowing for fast random access to any subsequence in the indexed FASTA file, 
while loading a minimal amount of the file in to memory. 
[pyfaidx](https://pypi.org/project/pyfaidx/) module implements pure Python classes for indexing, retrieval, 
and in-place modification of FASTA files using a samtools compatible index.

### 2. maxentpy
[maxentpy](https://github.com/kepbod/maxentpy) is a python wrapper for MaxEntScan to calculate splice site strength.
It contains two functions. score5 is adapt from [MaxEntScan::score5ss](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) to score 5' splice sites. 
score3 is adapt from [MaxEntScan::score3ss](http://genes.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq_acc.html) to score 3' splice sites. 

maxentpy is already included in the **autopvs1**.

### 3. pyhgvs
[pyhgvs](https://github.com/counsyl/hgvs) provides a simple Python API for parsing, formatting, and normalizing HGVS names.

But it only supports python2, I modified it to support python3 and added some other features. It is also included in 
the **autopvs1**.

### 4. configure file
```bash
# autopvs1/config.ini
[DEFAULT]

ref = data/hg19.fa
trans = data/refGenePrimary.gpe

morbidmap = data/morbidmap2.txt
domain = data/PM1.well-established.domain.bed
hotspot = data/PM1.hotspot.bed
curated_region = data/PM1.expert_curated.bed
pathogenic_ref = data/clinvar_20190624_pathogenic.vcf
pvs1levels = data/PVS1.level
ba1_exception = data/BA1.exception
exon_lof_frequent = data/exon.lof_frequent.bed

```

## USAGE
### single variant
```python
from autopvs1.autoPVS1 import AutoPVS1
test1=AutoPVS1('22-36678800-G-A')
print(test1.consequence, test1.hgvs_c, test1.hgvs_p, test1.strength_raw)
```

### bath processing
VCF input are required
```bash
cd test
bash vep.sh pilot/pvs1.pilot.vcf >pvs1.pilot.vep
python3 vep_lof_filter.py pvs1.pilot.vep
python3 test.py pvs1.pilot.vep.lof >pvs1.pilot.vep.lof.autopvs1

```

#Terms of use
Users may freely use the AutoPVS1 for non-commercial purposes as long as they properly cite it. 
This resource is intended for research purposes only. For clinical or medical use, please consult professionals.