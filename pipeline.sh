#specify threadnum
TNUM=12
REF=GCF_000005575.2_AgamP3_genomic.fna

for f in `ls raw_reads/*.fastq.gz | sed 's/_[^_]*$//' |sed 's/_[^_]*$//' | sed 's!.*/!!' | uniq`
do
mkdir results/"$f"
TRIMDIR=results/"$f"/trimmed_reads
mkdir $TRIMDIR


#run fastp to trim raw reads
fastp -w $TNUM -i raw_reads/$f"_R1_001.fastq.gz" -I raw_reads/$f"_R2_001.fastq.gz" -o results/$f/trimmed_reads/$f"_1_trim.fq.gz" -O results/$f/trimmed_reads/$f"_2_trim.fq.gz"
	
#run fastqc on both pairs
fastqc -t $TNUM -o results/$f -f fastq $TRIMDIR/$f"_2_trim.fq.gz"
fastqc -t $TNUM -o results/$f -f fastq $TRIMDIR/$f"_1_trim.fq.gz"

#bwa, fixmate, sort, markdup
mkdir results/"$f"/tmp
bwa mem -t $TNUM ref/$REF results/$f/trimmed_reads/$f"_1_trim.fq.gz" results/$f/trimmed_reads/$f"_2_trim.fq.gz" | \
    samtools fixmate -m - - | \
    samtools sort -@$TNUM -T results/"$f"/tmp - | \
    samtools markdup -@$TNUM --reference $REF - results/"$f"/$f".srt.dp.bam"

#flagstat
samtools flagstat results/"$f"/$f".srt.dp.bam" > results/"$f"/$f".flagfile"

#qualimap
qualimap bamqc -bam results/"$f"/$f".srt.dp.bam"

done

#multiqc report generate
multiqc .

#collect flagstats
for f in results/*/*.flagfile
   do
   flag=`< \$f cut -d \\+ -f 1 | tr -s '[:blank:]' ','`
   echo "\${f%.stats.txt}",\$flag | tr -d ' ' >> results/mapping_stats.csv
   done
