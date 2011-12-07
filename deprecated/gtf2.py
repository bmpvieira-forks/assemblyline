'''
Created on May 4, 2010

@author: mkiyer
'''

GTF_EMPTY_FIELD = '.'
GTF_ATTR_SEP = ';'
GTF_ATTR_TAGVALUE_SEP = ' '

class GTF2Feature(object):
    '''
    1. seqname - The name of the sequence. Must be a chromosome or scaffold.
    2. source - The program that generated this feature.
    3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    5. end - The ending position of the feature (inclusive).
    6. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
    7. strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
    8. frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.

    chr1    Cufflinks       transcript      3109    3184    1000    .       .       gene_id "3064YAAXX_2.1"; transcript_id "3064YAAXX_2.1.0"; RPKM "0.1802096261"; frac "1.000000"; conf_lo "0.000000"; conf_hi "0.000000"; cov "0.108248";
    '''
    __slots__ = ('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attrs', 
                 'gene_id', 'transcript_id') 

    def __str__(self):
        line = [self.seqid,
                self.source,
                self.type,
                str(self.start),
                str(self.end),
                str(self.score),
                str(self.strand),
                self.phase]
        attr_str = ' '.join('%s "%s";' % (k, v) for (k, v) in self.attrs.iteritems())
        line.append(attr_str)
        return '\t'.join(line)
        
    @staticmethod
    def from_str(line, attr_defs=None):        
        f = GTF2Feature()        
        # read the GTF line
        fields = line.strip().split('\t')
        f.seqid = fields[0]
        f.source = fields[1]
        f.type = fields[2]
        f.start = int(fields[3])
        f.end = int(fields[4])
        f.score = 0 if (fields[5] == '.') else float(fields[5])
        strand = fields[6]
        if not (strand == '+' or strand == '-'):
            strand = GTF_EMPTY_FIELD
        f.strand = strand        
        f.phase = fields[7]        
        attrs = {}
        if fields[8] != GTF_EMPTY_FIELD:
            attr_strings = fields[8].split(GTF_ATTR_SEP)
            for a in attr_strings:
                a = a.strip()
                if len(a) == 0:
                    continue
                tag, value = a.split(GTF_ATTR_TAGVALUE_SEP, 1)
                # remove quotes
                value = value.split('"')[1]
                # apply parsing function
                if (attr_defs != None) and (tag in attr_defs) and (attr_defs[tag] != None):
                    value = attr_defs[tag](value)
                attrs[tag] = value
        f.attrs = attrs
        return f

def read_gtf2_file(fileh, attr_defs=None):
    for line in fileh:
        # read the GTF line
        if not line:
            continue
        yield GTF2Feature.from_str(line, attr_defs)

def _make_bed_line(features):
    chrom = None
    estarts = []
    eends = []
    name = None
    strand = None
    block_starts = []
    block_sizes = []

    for f in features:
        if (chrom == None):
            chrom = f.seqid
            name = f.attrs['transcript_id']
            strand = '.' if f.strand == None else f.strand
        else:
            assert f.seqid == chrom
            assert f.attrs['transcript_id'] == name
        estarts.append(f.start)
        eends.append(f.end)

    tx_start = min(estarts)
    tx_end = max(eends)
    for estart, eend in zip(estarts, eends):
        block_sizes.append(eend - estart)
        block_starts.append(estart - tx_start)
    # construct BED feature
    s = [chrom, str(tx_start), str(tx_end), 
         name, '0', strand,
         str(tx_start), str(tx_end),
         '.',
         str(len(estarts))]
    s.append(','.join(map(str, block_sizes)) + ',')
    s.append(','.join(map(str, block_starts)) + ',')
    return '\t'.join(s)

    
def gtf2bed(gtf_fileh, attr_defs=None):
    cur_tx_id = None
    cur_tx_features = []

    for feature in read_gtf2_file(gtf_fileh, attr_defs=attr_defs):
        tx_id = feature.attrs['transcript_id']
        if cur_tx_id == None:
            cur_tx_id = tx_id
        elif cur_tx_id != tx_id:
            bedstr = _make_bed_line(cur_tx_features)
            yield bedstr
            # reset current transcript
            cur_tx_id = tx_id
            cur_tx_features = []
        if feature.type == 'exon':
            cur_tx_features.append(feature)
    
    if len(cur_tx_features) > 0:
        bedstr = _make_bed_line(cur_tx_features)
        yield bedstr

    
    