'''
Created on Nov 28, 2011

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import logging

class LibraryInfo(object):
    fields = ('cohort_id', 'patient_id', 'sample_id', 'library_id',
              'gtf_file', 'bam_file')

    def __init__(self):
        for f in LibraryInfo.fields:
            setattr(self, f, None)

    @staticmethod
    def from_fields(fields, field_dict=None):
        if field_dict is None:
            field_dict = dict((x,i) for i,x in enumerate(LibraryInfo.fields))
        tbl = LibraryInfo()
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
            yield LibraryInfo.from_fields(fields, field_dict)
        fh.close()
        
    def is_valid(self):
        if self.gtf_file is None:
            logging.error("Library %s GTF file is None" % (self.library_id))
            return False
        if not os.path.exists(self.gtf_file):
            logging.error("Library %s GTF file not found" % (self.library_id))
            return False
        if self.bam_file is None:
            logging.error("Library %s BAM file is None" % (self.library_id))
            return False
        if not os.path.exists(self.bam_file):
            logging.error("Library %s BAM file not found" % (self.library_id))
            return False
        return True
