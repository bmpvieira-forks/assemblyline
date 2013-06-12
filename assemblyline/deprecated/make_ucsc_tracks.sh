PREFIX=$1
CHROM_SIZES_FILE=$2
BASEURL=$3

python ~/workspace/assemblyline/assemblyline/utils/ucsc_to_ensembl.py -r ${PREFIX}.bed > ${PREFIX}.ucsc.bed
python ~/workspace/assemblyline/assemblyline/utils/ucsc_to_ensembl.py -r ${PREFIX}_none.bedgraph > ${PREFIX}_none.ucsc.bedgraph
python ~/workspace/assemblyline/assemblyline/utils/ucsc_to_ensembl.py -r ${PREFIX}_neg.bedgraph > ${PREFIX}_neg.ucsc.bedgraph
python ~/workspace/assemblyline/assemblyline/utils/ucsc_to_ensembl.py -r ${PREFIX}_pos.bedgraph > ${PREFIX}_pos.ucsc.bedgraph

# BED file
sed 1,1d ${PREFIX}.ucsc.bed > ${PREFIX}.ucsc.noheader.bed
sort -k1,1 -k2,2n ${PREFIX}.ucsc.noheader.bed > ${PREFIX}.ucsc.srt.bed
bedToBigBed ${PREFIX}.ucsc.srt.bed ${CHROM_SIZES_FILE} ${PREFIX}.ucsc.srt.bb
TRACK_LINE=`head -n 1 ${PREFIX}.ucsc.bed`
echo "$TRACK_LINE bigDataUrl=http://${BASEURL}/${PREFIX}.ucsc.srt.bb"
rm ${PREFIX}.ucsc.noheader.bed

# BEDGRAPH FILES
python ~/workspace/assemblyline/assemblyline/utils/bedgraph_to_bigwig.py --baseurl ${BASEURL} ${PREFIX}_none.ucsc.bedgraph ${CHROM_SIZES_FILE}
python ~/workspace/assemblyline/assemblyline/utils/bedgraph_to_bigwig.py --baseurl ${BASEURL} ${PREFIX}_pos.ucsc.bedgraph ${CHROM_SIZES_FILE}
python ~/workspace/assemblyline/assemblyline/utils/bedgraph_to_bigwig.py --baseurl ${BASEURL} ${PREFIX}_neg.ucsc.bedgraph ${CHROM_SIZES_FILE}

