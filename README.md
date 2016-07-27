# OVarCall script for python package

# Motivation
Detection of somatic SNVs of low allele frequency (under 7%) in exome sequence data is stil difficult.
Overlapping paired-end read is reported to be important for detection of mutations of low allele frequency in PCR targeted deep sequence, so usage of these information is expected to be important.
However, existing methods for PCR targeted deep sequence data use only overlapping paired-end reads and was not designed for usual exome sequence data that include overlapping and non-overlapping paired-end reads.
We constructed a Bayesian hierarchical method, OVarCall, for the detection of somatic mutations with low allele frequencies from exome sequence data.


# Paper
OVarCall: Bayesian mutation calling method utilizing overlapping paired-end reads, Proc. the 12th International Symposium on Bioinformatics Research and Applications, Lecture Notes in Bioinformatics, Springer-Verlag Berlin Heigelberg, 9683, 40-51, 2016.


# Dependency
## Software
samtools 

## Python
Python(>=2.7), pysam, scipy,pyVCF


##  Install
```
git clone https://github.com/takumorizo/OVarCall
cd OVarCall
python setup.py build
python setup.py install
```


## Run
You can use OVarCall from Tumor bam and Normal bam.
```
OVarCall [-h] [--version] -1 BAM1 -2 BAM2 -o OUTPUT -r REF_FA -s SAMTOOLS_PATH -p PARAMETER_SETTINGS [-R REGION] [-l LOG_LEVEL]

-1: Input tumor bam path.
-2: Input normal bam path.
-o: Output file path.
-r: Reference fasta file path.
-s: Samtools path.
-p: Parameter settings path , ex) ./OVarCall.ini.
-R: Genomic regions , ex) chr1:1000-2000.
-l: Loglevel, ex) CRITICAL,ERROR,WARNING,INFO,DEBUG.
```

You can also use OVarCall from analysed vcf file or annovar file, and scores are added.
```
OVarFilter [-h] [--version] -1 BAM1 -2 BAM2 -i INPUT -o OUTPUT -r REF_FA -p PARAMETER_SETTINGS [-f] [-l LOG_LEVEL] [-a]

-1: Input tumor bam path.
-2: Input normal bam path.
-i; Input vcf or annovar file path.
-a: set this option when you set annovar file in -i.
-o: Output file path.
-f: Activate pileup filter.
-r: Reference fasta file path.
-p: Parameter settings path , ex) ./OVarCall.ini.
-l: Loglevel, ex) CRITICAL,ERROR,WARNING,INFO,DEBUG.
```

## Output format
The meaning of output is as follows.
```
TYPE : M(SNV), I(Insertion), D(Deletion)
Chr : Chromosome
pos : 1-indexed position
ref : base @ position
obs : altered base(TYPE == M), inserted base(TYPE == I), deleted base(TYPE == D).
scoreWithOverlap : log_10_Bayes_factor, considering overlap
scoreWithoutOverlap : log_10_Bayes_factor, ignoring overlap
TumorPileup, NormalPileup: Ref+,Obs+,Other+,Ref-,Obs-,Other-,RefRef,RefObs,RefOther,ObsRef,ObsObs,ObsOther,OtherRef,OtherObs,OtherOther
```