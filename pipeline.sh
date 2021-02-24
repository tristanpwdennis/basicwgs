#init conda
source ~/.bashrc
cd /export/projects/III-data/lamberton/tdd3v/marshallagia/basicwgs

#specify threadnum
TNUM=12
REF=tc_canu1.9_200521.fa

rm -rf results
mkdir results

for f in `ls raw_reads/*.fq.gz | sed 's/_[^_]*$//' | sed 's!.*/!!' | uniq`
do
mkdir results/"$f"
TRIMDIR=results/"$f"/trimmed_reads
mkdir $TRIMDIR


#run fastp to trim raw reads
fastp -w $TNUM -i raw_reads/$f"_1.fq.gz" -I raw_reads/$f"_2.fq.gz" -o results/$f/trimmed_reads/$f"_1_trim.fq.gz" -O results/$f/trimmed_reads/$f"_2_trim.fq.gz"
	
#run fastqc on both pairs
fastqc -t $TNUM -o results/$f -f fastq $TRIMDIR/$f"_2_trim.fq.gz"
fastqc -t $TNUM -o results/$f -f fastq $TRIMDIR/$f"_1_trim.fq.gz"

#bwa, fixmate, sort, markdup
mkdir results/"$f"/tmp
bwa mem -t $TNUM ref/$REF results/$f/trimmed_reads/$f"_1_trim.fq.gz" results/$f/trimmed_reads/$f"_2_trim.fq.gz" | \
    samtools fixmate -u -m - - | \
    samtools sort -u -@$TNUM -T results/"$f"/tmp - | \
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
