'''
Created on Dec 10, 2012

@author: mkiyer
'''

def get_total_reads(alignment_metrics_file):
    f = open(alignment_metrics_file)
    for line in f:
        if line.startswith("## METRICS"):
            break
    header_fields = f.next().strip().split('\t')
    total_reads = 0
    for line in f:
        fields = line.strip().split('\t')
        if fields[0] == "FIRST_OF_PAIR":
            total_reads += int(fields[header_fields.index("TOTAL_READS")])
        elif fields[0] == "SECOND_OF_PAIR":
            total_reads += int(fields[header_fields.index("TOTAL_READS")])
    return total_reads
