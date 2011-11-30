'''
Created on Apr 19, 2010

@author: mkiyer
'''
from bx.intervals.cluster import ClusterTree

class Transcript(object):
    __slots__ = ('id', 'chrom', 'strand', 'score', 'exons')
    
    def __init__(self, id=None, chrom=None, strand=None, score=0.0, exons=None):
        self.id = id
        self.chrom = chrom
        self.strand = '.' if strand == None else strand
        self.score = score
        self.exons = [] if exons == None else exons

    def __str__(self):
        return ("<Transcript id=%s chrom=%s strand=%s score=%.3f exons=%s" %
                (self.id, self.chrom, self.strand, self.score, self.exons))

    @property
    def tx_start(self):
        return self.exons[0][0]
    @property
    def tx_end(self):
        return self.exons[-1][1]
    @property
    def introns(self):
        return zip([e[1] for e in self.exons[:-1]],
                   [e[0] for e in self.exons[1:]])

    @staticmethod
    def from_gtf_features(features, score_attr=None, exon_score=0):
        score = 0.0
        if score_attr != None:            
            score = sum(f.attrs[score_attr] for f in features)
        if exon_score > 0:
            score += exon_score * len(features)
        # cufflinks has some strandedness bugs in single-exon transcripts
        # so never set the strand of 1-exon transcripts
        if len(features) == 1:
            strand = '.'
        else:
            strand = features[0].strand
        t = Transcript(id=features[0].attrs['transcript_id'],
                       chrom=features[0].seqid,
                       strand=strand,
                       score=score,
                       exons=[(f.start, f.end) for f in features])
        return t

    @staticmethod
    def merge(tx1, tx2, id):
        mergetx = Transcript(id, tx1.chrom, tx1.strand)
        assert tx1.chrom == tx2.chrom
        assert tx1.strand == tx2.strand
        # need to cluster the exons to handle
        tree = ClusterTree(0, 1)
        for e in tx1.exons:
            tree.insert(e[0], e[1], tx1.id)
        for e in tx2.exons:
            tree.insert(e[0], e[1], tx2.id)
        for cluster in tree.getregions():
            start, end, ids = cluster
            mergetx.exons.append((start,end))
        return mergetx

    def to_bed12(self):
        block_sizes = []
        block_starts = []
        for e0, e1 in self.exons:
            block_starts.append(e0 - self.tx_start)
            block_sizes.append(e1 - e0)        
        # write
        s = '\t'.join([self.chrom, 
                       str(self.tx_start), 
                       str(self.tx_end),
                       str(self.id),
                       str(self.score),
                       self.strand if self.strand != None else '.',
                       str(self.tx_start),
                       str(self.tx_start),
                       '0',
                       str(len(self.exons)),
                       ','.join(map(str,block_sizes)) + ',',
                       ','.join(map(str,block_starts)) + ','])
        return s
