hello Amplicon!


# WORK FLOW


## Reads count

```
 files=`ls ./fastq/*_S*_L001_R1_001.fastq.gz`
 ${HOME}/biotools/local/assembly/bin/seqkit stats $files > out.txt

 perl -F'\s+' -anle 'BEGIN{ $info={}; } next if (/^file/);  $file= $1 if ($F[0]=~/.*\/(\S+)$/); $lib= $1 if ($file=~/^(\S+)-P_S\d+_L001_R1_001.fastq.gz/); print "$lib\t$file\t$F[3]\t$F[4]\t$F[6]"; END{ }' out.txt

```
