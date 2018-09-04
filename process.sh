y=$1
BINDIR=$2
OUTDIR=$3
bamFile=$4
fastaFile=$5
echo 'flagsflags: ' $1 ' ' $2 ' ' $3 ' ' $4 ' ' $5
numReads=`wc -l < $y`
l=${#y}
ol=${#OUTDIR}
ty=${y:ol+9:l}
c="${ty#*.}"
c="${c#*.}"
x=${y:ol+9:l-4-ol-9-${#c}-1}
echo 'Processing insertion '$c:$x

# Produce an empty file to indicate this insertion has started being processed
touch $OUTDIR/cert/$x'.txt.'$c'.cert'

# Produce fastq file of reads supporting the insertion based on alignments
samtools view -H $bamFile > $OUTDIR/deletes/tmp.txt
samtools view $bamFile  $c:$(($x - 10000))-$(($x + 10000)) | grep -w -f $OUTDIR/deletes/$x.txt.$c > $OUTDIR/deletes/supportingreadalignments_$x.sam.$c
cat $OUTDIR/deletes/supportingreadalignments_$x.sam.$c >> $OUTDIR/deletes/tmp.txt
mv $OUTDIR/deletes/tmp.txt $OUTDIR/deletes/supportingreadalignments_$x.sam.$c
java -cp "${BINDIR}" ReadInsertExtraction $OUTDIR/deletes/supportingreadalignments_$x.sam.$c $x > $OUTDIR/deletes/supportingreadswithinsertions_$x.fa.$c
samtools bam2fq $OUTDIR/deletes/supportingreadalignments_$x.sam.$c > $OUTDIR/deletes/supportingreads_$x.fq.$c

# Convert from fastq to fasta format
java -cp "${BINDIR}" FastaFileFixer $OUTDIR/deletes/supportingreads_$x.fq.$c $OUTDIR/deletes/supportingreads_$x.fa.$c
rm $OUTDIR/deletes/fixed_$x.fq.$c

# Put the reads in a format readable by FalconSense
java -cp "${BINDIR}" FalconFormatter $OUTDIR/deletes/supportingreads_$x.fa.$c $OUTDIR/deletes/supportingreads_falconsenseinput_$x.fa.$c
java -cp "${BINDIR}" FalconFormatter $OUTDIR/deletes/supportingreadswithinsertions_$x.fa.$c $OUTDIR/deletes/supportingreadswithinsertions_falconsenseinput_$x.fa.$c

# Run Falconsense using default parameters
quarterReads=`expr $numReads / 8`
$BINDIR/falcon_sense --min_idt 0.7 --min_len 500 --max_read_len 1234567 --min_ovl_len 250 --min_cov $quarterReads --n_core 2 < $OUTDIR/deletes/supportingreads_falconsenseinput_$x.fa.$c > $OUTDIR/falconsense_output/$x.correctedReads.fasta.$c
$BINDIR/falcon_sense --min_idt 0.7 --min_len 500 --max_read_len 1234567 --min_ovl_len 250 --min_cov 2 --n_core 2 < $OUTDIR/deletes/supportingreadswithinsertions_falconsenseinput_$x.fa.$c > $OUTDIR/falconsense_output/$x.insertions_correctedReads.fasta.$c

# Reformat the file to have the entire sequences on one line for downstream processing
java -cp "${BINDIR}" FastaFileFixer2 $OUTDIR/falconsense_output/$x.correctedReads.fasta.$c $OUTDIR/falconsense_output/$x.correctedReads.fixed.fasta.$c
mv $OUTDIR/falconsense_output/$x.correctedReads.fixed.fasta.$c $OUTDIR/falconsense_output/$x.correctedReads.fasta.$c
java -cp "${BINDIR}" FastaFileFixer2 $OUTDIR/falconsense_output/$x.withinsertions_correctedReads.fasta.$c $OUTDIR/falconsense_output/$x.insertions_correctedReads.fixed.fasta.$c
mv $OUTDIR/falconsense_output/$x.insertions_correctedReads.fixed.fasta.$c $OUTDIR/falconsense_output/$x.insertions_correctedReads.fasta.$c

offset=100000
left=$(($x - $offset))
if [ "$left" -lt 1 ]
then
    left=1
fi
echo 'Aligning assembled supporting reads back to: ' $c:$left-$(($x + $offset))
echo 'Reference genome file: ' $fastaFile
samtools faidx $fastaFile $c:$left-$(($x + $offset)) > $OUTDIR/reference_snippets/"$x".fa.$c

ngmlr -t 4 -r $OUTDIR/reference_snippets/"$x".fa.$c -q $OUTDIR/falconsense_output/"$x".correctedReads.fasta.$c -o $OUTDIR/alignment_of_assembled_reads/"$x".sam.$c
ngmlr -t 4 -r $OUTDIR/reference_snippets/"$x".fa.$c -q $OUTDIR/falconsense_output/"$x".insertions_correctedReads.fasta.$c -o $OUTDIR/alignment_of_assembled_reads/insertions_"$x".sam.$c
echo 'aligned reads'
delete=`java -cp "${BINDIR}" BestDeletionFinder $OUTDIR/alignment_of_assembled_reads/"$x".sam.$c $x $offset 'LEN'`
pos=`java -cp "${BINDIR}" BestDeletionFinder $OUTDIR/alignment_of_assembled_reads/"$x".sam.$c $x $offset 'POS'`

delete2=`java -cp "${BINDIR}" BestDeletionFinder $OUTDIR/alignment_of_assembled_reads/insertions_"$x".sam.$c $x $offset 'LEN'`
pos2=`java -cp "${BINDIR}" BestDeletionFinder $OUTDIR/alignment_of_assembled_reads/insertions_"$x".sam.$c $x $offset 'POS'`

echo $x >> $OUTDIR/refined_deletion_list.txt
echo 'delete: '$delete >> $OUTDIR/refined_deletion_list.txt

echo '>delete_'$c':'$x > $OUTDIR/seqs/$x.$c.fa
echo '>delete_'$c':'$x > $OUTDIR/seqs/$x.$c.pos
if [ "$delete" != '' ]
then
    echo $delete >> $OUTDIR/seqs/$x.$c.fa
    echo $pos >> $OUTDIR/seqs/$x.$c.pos
else
    echo $delete2 >> $OUTDIR/seqs/$x.$c.fa
    echo $pos2 >> $OUTDIR/seqs/$x.$c.pos
fi

rm $OUTDIR/reference_snippets/"$x".*

# Produce a file to indicate this deletion has finished being processed
touch $OUTDIR/cert/$x'.txt.'$c'.cert.done'
