'''
Created on Nov 22, 2011

@author: mkiyer
'''

def parse_raw_juncs(line_iter):
    """
    generator to parse lines in the tophat "raw juncs" format
    """
    for line in line_iter:
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[3]
        yield (chrom, start, end, strand)
        
