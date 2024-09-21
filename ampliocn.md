hello Amplicon!


# WORK FLOW

## Fastq format

![Fastq image](images/fastq_code.png)



## Reads count

> [!Tip]
> Count of reads
```
 files=`ls ./fastq/*_S*_L001_R1_001.fastq.gz`
 ${HOME}/biotools/local/assembly/bin/seqkit stats $files > out.txt

 perl -F'\s+' -anle 'BEGIN{ $info={}; } next if (/^file/);  $file= $1 if ($F[0]=~/.*\/(\S+)$/); $lib= $1 if ($file=~/^(\S+)-P_S\d+_L001_R1_001.fastq.gz/); print "$lib\t$file\t$F[3]\t$F[4]\t$F[6]"; END{ }' out.txt

```


1. COMBINING PAIRED-END

> [!Note]
> RAWデータが"fastqディレクトリー"にある前提

```
 mkdir pear; \
  for i in ./fastq/*_L001_R1_001.fastq.gz; do a=$i; \
   j=${i##./*/}; k=${j%%-P_S*_L001_R1_001.fastq.gz}; echo "$k";
  $HOME/biotools/local/assembly/pear-0.9.10-bin-64/pear-0.9.10-bin-64 \
    -f $i -r ${i/_R1/_R2} -o ./pear/${k}_pear -j 12; \
 done

```
OUTPUT: ./pear/*_pear.assembled.fastq

2. REMOVE PRIMER
    -V4-V5

```
export PATH="/home/impact/biotools/rhel6/miniconda3/bin:$PATH"
mkdir cleanup; \
 for i in ./pear/*_pear.assembled.fastq; do a=$i; \
    j=${i##./*/}; k=${j%%_pear.assembled.fastq}; a=$i; echo "$k"; \
 cutadapt -j 12 -e 0.1 -g file:$HOME/Desktop/work_pop/bin/primer_fwd_iupac.fas -n 5 \
   -o out.fastq $i; \
 cutadapt -j 12 -e 0.1 -a file:$HOME/Desktop/work_pop/bin/primer_rev_iupac_comp.fas -n 5 \
   -o ./cleanup/${k}_pear_noprim.fq out.fastq; \
 done
```

OUTPUT: ./cleanup/*_pear_noprim.fq


3. FILTER QC & LENGTH
   

```
for i in ./cleanup/*_pear_noprim.fq; do a=$i; \
    j=${i##./*/}; k=${j%%_pear_noprim.fq}; a=$i; echo "$k"; \
 $HOME/Desktop/work_pop/bin/pickup_qc_fastq.pl -se ./cleanup/${k}_pear_noprim.fq \
        -qc 30 0.98 -minlen 100 -maxlen 550; \
 mv out_pickup_R1_SE.fq ./cleanup/${k}_pear_noprim_qc.fq; \
 done

```

OUTPUT: ./cleanup/*_pear_noprim_qc.fq


> [!Tip]
> Count of reads

```
 files=`ls ./cleanup/*_pear_noprim_qc.fq`
 ${HOME}/biotools/local/assembly/bin/seqkit stats $files > out.txt

 perl -F'\s+' -anle 'BEGIN{ $info={}; } next if (/^file/);  $file= $1 if ($F[0]=~/.*\/(\S+)$/); $lib= $1 if ($file=~/^(\S+)_pear_noprim_qc.fq/); print "$lib\t$file\t$F[3]\t$F[4]\t$F[6]"; END{ }' out.txt

```
