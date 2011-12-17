'''
Created on Nov 28, 2011

@author: mkiyer
'''

class LaneInfo(object):
    __slots__ = ('cohort', 'patient', 'sample', 'library', 'lane', 'qc',
                 'use_juncs', 'use_transcripts',  'valid', 
                 'aligned_reads', 'read_length', 'tophat_juncs_file',
                 'cufflinks_gtf_file')

    @staticmethod
    def from_fields(fields, field_dict=None):
        if field_dict is None:
            field_dict = dict((x,i) for i,x in enumerate(LaneInfo.__slots__))
        l = LaneInfo()
        for attrname in LaneInfo.__slots__:
            setattr(l, attrname, fields[field_dict[attrname]])
        # convert types        
        l.use_juncs = int(l.use_juncs)
        l.use_transcripts = int(l.use_transcripts)
        l.valid = int(l.valid)
        if l.aligned_reads == "NA":
            l.aligned_reads = None
        else:
            l.aligned_reads = int(float(l.aligned_reads))
        l.read_length = int(l.read_length)
        if l.tophat_juncs_file == "NA": 
            l.tophat_juncs_file = None
        if l.cufflinks_gtf_file == "NA": 
            l.cufflinks_gtf_file = None
        return l

    def is_valid(self):
        if self.qc == "FAIL":
            return False
        if self.aligned_reads is None:
            return False
        if self.use_juncs and (self.tophat_juncs_file is None):
            return False
        if self.use_transcripts and (self.cufflinks_gtf_file is None):
            return False
        return True

    @staticmethod
    def from_file(filename):
        fh = open(filename)
        # header
        field_names = fh.next().strip().split('\t')
        field_dict = dict((x,i) for i,x in enumerate(field_names))
        # table rows
        for line in fh:
            fields = line.strip().split('\t')
            yield LaneInfo.from_fields(fields, field_dict)
        fh.close()
 
class LibInfo(object):
    __slots__ = ('cohort', 'patient', 'sample', 'library', 'lanes',
                 'library_type', 'frag_len_mean', 'frag_len_std_dev',
                 'has_pe_lanes', 'valid', 'for_assembly', 
                 'cufflinks_gtf_file', 'bam_file') 

    @staticmethod
    def from_fields(fields, field_dict=None):
        if field_dict is None:
            field_dict = dict((x,i) for i,x in enumerate(LibInfo.__slots__))
        l = LibInfo()
        for attrname in LibInfo.__slots__:
            setattr(l, attrname, fields[field_dict[attrname]])
        # convert types        
        if l.frag_len_mean == "NA":
            l.frag_len_mean = None
        else:
            l.frag_len_mean = int(l.frag_len_mean)
        if l.frag_len_std_dev == "NA":
            l.frag_len_std_dev = None
        else:
            l.frag_len_std_dev = int(l.frag_len_std_dev)
        l.valid = int(l.valid)
        l.for_assembly = int(l.for_assembly)
        l.has_pe_lanes = int(l.has_pe_lanes)
        if l.cufflinks_gtf_file == "NA": 
            l.cufflinks_gtf_file = None
        if l.bam_file == "NA": 
            l.bam_file = None
        return l

    def is_valid(self):
        if not self.valid:
            return False
        if self.for_assembly and (self.cufflinks_gtf_file is None):
            return False
        if self.bam_file is None:
            return False
        return True

    @staticmethod
    def from_file(filename):
        fh = open(filename)
        # header
        field_names = fh.next().strip().split('\t')
        field_dict = dict((x,i) for i,x in enumerate(field_names))
        # table rows
        for line in fh:
            fields = line.strip().split('\t')
            yield LibInfo.from_fields(fields, field_dict)
        fh.close()
