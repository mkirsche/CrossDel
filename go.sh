BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
number=$RANDOM ##makes it work for multiple instances.
mkdir $WORKINGDIR/csout_$number
OUTDIR=$WORKINGDIR/csout_$number
echo $WORKINGDIR
echo $BINDIR
echo $OUTDIR
rm -r $OUTDIR/*
rm -r $OUTDIR
export _JAVA_OPTIONS="-XX:ParallelGCThreads=2"
javac $BINDIR/*.java
mkdir $OUTDIR
mkdir $OUTDIR/falconsense_output
mkdir $OUTDIR/deletes
mkdir $OUTDIR/alignment_of_assembled_reads
mkdir $OUTDIR/seqs
mkdir $OUTDIR/reference_snippets
mkdir $OUTDIR/cert

usage() { echo "Usage: $0 -v <vcfFile> -b <bamFile> -f <fastaFile> -o <outputFile>" 1>&2; exit 1; }

while getopts v:b:f:o: option
do
    case "${option}"
    in
    v) vcfFile=${OPTARG};;
    b) bamFile=${OPTARG};;
    f) fastaFile=${OPTARG};;
    o) outputFile=$OPTARG;;
 esac
done

if [ -z "${vcfFile}" ] || [ -z "${bamFile}" ] || [ -z "${fastaFile}" ] || [ -z "${outputFile}" ]; then
    usage
fi

if [ ! -r $bamFile'.bai' ]
then
  echo "Indexing bam file"
  samtools index $bamFile
fi

vcfPath=$WORKINGDIR/$vcfFile
if [[ $vcfFile == /* ]]; then
    vcfPath=$vcfFile
fi

fastaPath=$WORKINGDIR/$fastaFile
if [[ $fastaFile == /* ]]; then
    fastaPath=$fastaFile
fi

echo 'Reference genome file: ' $fastaPath
echo 'VCF file to be refined: ' $vcfPath

DELETE_BEFORE=1 # The number of characters before the insertion to include in the REF field of the new VCF file
DELETE_AFTER=0 # The number of characters after the insertion to include in the REF field of the new VCF file

# Generate lists of reads for all insertions
java -cp "${BINDIR}" ReadFinder $vcfPath $OUTDIR/deletes

numSupportedVariants=`cat $OUTDIR/inserts/out.log`

if [ "$numSupportedVariants" = "0" ]; then
   echo "No variant with supporting reads found";
   exit;
fi

echo 'bin dir: '$BINDIR
echo 'out dir: '$OUTDIR
echo 'bam file: '$bamFile
echo 'fasta path: '$fastaPath
numFiles=`ls $OUTDIR/deletes/*.txt.* | wc -l`
echo 'number of deletions to process: '$numFiles
# Process all insertions in parallel
parallel --gnu --timeout 500 --jobs 16 "${BINDIR}"/process.sh {} $BINDIR $OUTDIR $bamFile $fastaPath ::: $OUTDIR/deletes/*.txt.*

wait

"${BINDIR}"/clean_parallel.sh $BINDIR $OUTDIR $bamFile $fastaFile
cat $OUTDIR/seqs/*.fa > $OUTDIR/all.seq
cat $OUTDIR/seqs/*.pos > $OUTDIR/all.pos
java -cp "${BINDIR}" VCFEditor $OUTDIR/all.seq $OUTDIR/all.pos $vcfPath $fastaPath $WORKINGDIR/$outputFile $DELETE_BEFORE $DELETE_AFTER


