'''
Created on Dec 10, 2012

@author: mkiyer
'''

class AlignmentSummaryMetrics(object):
    fields = ["CATEGORY",
              "TOTAL_READS",
              "PF_READS",
              "PCT_PF_READS",
              "PF_NOISE_READS",
              "PF_READS_ALIGNED",
              "PCT_PF_READS_ALIGNED",
              "PF_ALIGNED_BASES",
              "PF_HQ_ALIGNED_READS",
              "PF_HQ_ALIGNED_BASES",
              "PF_HQ_ALIGNED_Q20_BASES",
              "PF_HQ_MEDIAN_MISMATCHES",
              "PF_MISMATCH_RATE",
              "PF_HQ_ERROR_RATE",
              "PF_INDEL_RATE",
              "MEAN_READ_LENGTH",
              "READS_ALIGNED_IN_PAIRS",
              "PCT_READS_ALIGNED_IN_PAIRS",
              "BAD_CYCLES",
              "STRAND_BALANCE",
              "PCT_CHIMERAS",
              "PCT_ADAPTER"]
    
    def __init__(self, alignment_metrics_file):
        f = open(alignment_metrics_file)
        for line in f:
            if line.startswith("## METRICS"):
                break
        self.category_field_dict = {}
        for line in f:
            fields = line.strip().split('\t')
            self.category_field_dict[fields[0]] = fields

    def get_total_reads(self):
        total_reads = 0
        index = AlignmentSummaryMetrics.fields.index("TOTAL_READS")
        for category,fields in self.category_field_dict.iteritems():
            if ((category == "FIRST_OF_PAIR") or 
                (category == "SECOND_OF_PAIR")):
                total_reads += int(fields[index])
        return total_reads

    def get_mean_read_length(self):
        index = AlignmentSummaryMetrics.fields.index("MEAN_READ_LENGTH")
        for category,fields in self.category_field_dict.iteritems():
            if ((category == "FIRST_OF_PAIR") or 
                (category == "SECOND_OF_PAIR")):
                return int(round(float(fields[index])))
        return None

def get_rnaseq_metrics(filename):
    f = open(filename)
    for line in f:
        if line.startswith("## METRICS"):
            fields = line.strip().split('\t')
            if fields[1] == "net.sf.picard.analysis.RnaSeqMetrics":
                break
    keys = f.next().strip().split('\t')
    vals = f.next().strip().split('\t')
    f.close()
    metrics_dict = dict(zip(keys,vals))
    return metrics_dict
    
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
