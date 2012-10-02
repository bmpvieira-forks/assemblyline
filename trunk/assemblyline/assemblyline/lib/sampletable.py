'''
Created on Nov 28, 2011

@author: mkiyer
'''
import os

class SampleInfo(object):
    fields = ('cohort', 'patient', 'sample', 'library',
              'gtf_file', 'bam_file')

    def __init__(self):
        self.cohort = None
        self.patient = None
        self.sample = None
        self.library = None
        self.gtf_file = None
        self.bam_file = None

    @staticmethod
    def from_fields(fields, field_dict=None):
        if field_dict is None:
            field_dict = dict((x,i) for i,x in enumerate(SampleInfo.fields))
        tbl = SampleInfo()
        for attrname in field_dict:
            setattr(tbl, attrname, fields[field_dict[attrname]])
        return tbl

    @staticmethod
    def from_file(filename):
        fh = open(filename)
        # header
        field_names = fh.next().strip().split('\t')
        field_dict = dict((x,i) for i,x in enumerate(field_names))
        # table rows
        for line in fh:
            fields = line.strip().split('\t')
            yield SampleInfo.from_fields(fields, field_dict)
        fh.close()
        
    def is_valid(self):
        if self.gtf_file is None:
            return False
        if not os.path.exists(self.gtf_file):
            return False
        if self.bam_file is None:
            return False
        if not os.path.exists(self.bam_file):
            return False
        return True
