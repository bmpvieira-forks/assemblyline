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
import xlrd
import xml.etree.cElementTree as etree
import collections

from assemblyline.rnaseq.lib.base import indent_xml

# fragment layouts
FRAGMENT_LAYOUT_SINGLE = "single"
FRAGMENT_LAYOUT_PAIRED = "paired"
FRAGMENT_LAYOUT_UNKNOWN = "unknown"
FRAGMENT_LAYOUTS = (FRAGMENT_LAYOUT_SINGLE, FRAGMENT_LAYOUT_PAIRED)

# strand protocols
FR_FIRSTSTRAND = 'fr-firststrand'
FR_UNSTRANDED = 'fr-unstranded'
LIBRARY_TYPES = (FR_FIRSTSTRAND, FR_UNSTRANDED)

# library protocols
EXOME_CAPTURE_TRANSCRIPTOME = 'exome_capture_transcriptome'
POLYA_TRANSCRIPTOME = 'polya_transcriptome'
LIBRARY_PROTOCOLS = (EXOME_CAPTURE_TRANSCRIPTOME,
                     POLYA_TRANSCRIPTOME)
# xml
LIBRARY_ROOT_TAG = "libraries"

class Library(object):
    fields = ('study_id', 'cohort_id', 'patient_id', 'sample_id', 
              'library_id', 'description', 'species', 'library_protocol', 
              'library_type', 'fragment_layout', 'seq_repo', 'read1_files', 
              'read2_files', 'bam_files')  

    def __init__(self, **kwargs):
        # basic parameters
        for attrname in Library.fields:
            setattr(self, attrname, kwargs.get(attrname, ''))
        # custom parameters
        self.params = kwargs["params"]
        # modify specific fields
        if not self.read1_files:
            self.read1_files = []
        else:
            self.read1_files = kwargs['read1_files'].split(",")              
        if not self.read2_files:
            self.read2_files = []
        else:        
            self.read2_files = kwargs['read2_files'].split(",")
        if not self.bam_files:
            self.bam_files = []
        else:
            self.bam_files = kwargs['bam_files'].split(",")
        # determine fragment layout
        if self.fragment_layout not in FRAGMENT_LAYOUTS:        
            if (len(self.read1_files) > 0 and
                len(self.read2_files) > 0):
                self.fragment_layout = FRAGMENT_LAYOUT_PAIRED
            elif (len(self.read1_files) > 0 and
                  len(self.read2_files) == 0):
                self.fragment_layout = FRAGMENT_LAYOUT_SINGLE
            else:
                self.fragment_layout = FRAGMENT_LAYOUT_UNKNOWN

    @staticmethod
    def from_fields(fields, field_dict=None):
        if field_dict is None:
            field_dict = dict((x,i) for i,x in enumerate(Library.fields))
        kwargs = {"params": {}}
        for attrname,ind in field_dict.iteritems():
            kwargs[attrname] = fields[ind].strip()
        return Library(**kwargs)

    @staticmethod
    def from_file(filename):
        fh = open(filename)
        # header
        field_names = fh.next().strip().split('\t')
        field_dict = dict((x,i) for i,x in enumerate(field_names))
        # table rows
        for line in fh:
            fields = line.strip().split('\t')
            yield Library.from_fields(fields, field_dict)
        fh.close()

    @staticmethod
    def from_xml_elem(elem):
        kwargs = {}
        for f in Library.fields:
            kwargs[f] = elem.findtext(f)
        params = {}
        for param_elem in elem.findall("param"):
            params[param_elem.get("name")] = param_elem.text
        kwargs["params"] = params
        return Library(**kwargs)

    @staticmethod
    def from_xml_file(xmlfile):
        tree = etree.parse(xmlfile)        
        root = tree.getroot()
        if root.tag == "library":
            yield Library.from_xml_elem(root)
        else:
            for elem in root.findall("library"):
                yield Library.from_xml_elem(elem)
    
    def to_xml(self, parent):
        root = etree.SubElement(parent, "library")
        remaining_fields = list(Library.fields)
        multi_fields = 'read1_files', 'read2_files', 'bam_files'
        for f in multi_fields:
            val = ','.join(getattr(self,f))
            elem = etree.SubElement(root, f)
            elem.text = val
            remaining_fields.remove(f)            
        for f in remaining_fields:
            elem = etree.SubElement(root, f)
            elem.text = getattr(self,f)
        for k,v in self.params.iteritems():
            elem = etree.SubElement(root, "param", name=k)
            elem.text = v
        return parent

    def to_xml_file(self, xmlfile):
        root = etree.Element(LIBRARY_ROOT_TAG)
        self.to_xml(root)
        # indent for pretty printing
        indent_xml(root)
        # write to file
        f = open(xmlfile, "w")
        print >>f, etree.tostring(root)
        f.close() 

    def is_valid(self):
        is_valid = True
        # check that RNA samples have valid strand protocol
        if (self.library_type not in LIBRARY_TYPES):
            logging.error("Invalid library type %s" % (self.library_type))
            is_valid = False
        if self.fragment_layout not in FRAGMENT_LAYOUTS:
            logging.error("Unknown fragment layout")
            is_valid = False    
        for filename in self.read1_files:
            if not os.path.exists(filename):
                logging.error("Library %s read 1 file %s not found" % (self.library_id, filename))
                is_valid = False
        if (self.fragment_layout == FRAGMENT_LAYOUT_PAIRED):
            if len(self.read1_files) != len(self.read2_files):
                logging.error("Library %s unequal number of read1 and read2 files" % (self.library_id))
            for filename in self.read2_files:
                if not os.path.exists(filename):
                    logging.error("Library %s read 2 file %s not found" % (self.library_id, filename))
                    is_valid = False
            for filename1 in self.read1_files:
                for filename2 in self.read2_files:
                    if filename1 == filename2:
                        logging.error("Library %s read 1 and read 2 point to same file %s" % (self.library_id, filename1))
                        is_valid = False
        for filename in self.bam_files:
            if not os.path.exists(filename):
                logging.error("Library %s bam file %s not found" % (self.library_id, filename))
                is_valid = False
        return is_valid

def read_wksheet(wksheet):
    field_names = wksheet.row_values(0)
    for rownum in xrange(1, wksheet.nrows):
        fields = wksheet.row_values(rownum)
        #print fields
        fields = [' '.join(str(field).split('\n')) for field in fields]
        # build dictionary of field names to field values
        field_name_value_dict = dict((field_names[i], fields[i]) for i in xrange(len(fields)))
        yield field_name_value_dict

def read_library_table_xls(filename):
    """
    parses an XLS file and constructs 'Patient', 'Sample', 'Library', 
    and 'Lane' objects stored as dictionaries keyed by unique id 
    """
    if not os.path.isfile(filename):
        raise OSError("File %s not found or not a regular file" % (filename))
    wkbook = xlrd.open_workbook(filename)
    # check that required sheet names exist
    sheet_names = wkbook.sheet_names()
    if not LIBRARY_ROOT_TAG in sheet_names:
        raise Exception("XLS file missing 'libraries' Sheet")
    if not "parameters" in sheet_names:
        raise Exception("XLS file missing 'parameters' Sheet")
    # read parameters
    params = {"patient": collections.defaultdict(lambda: {}),
              "sample": collections.defaultdict(lambda: {}),
              "library": collections.defaultdict(lambda: {})}
    for field_dict in read_wksheet(wkbook.sheet_by_name("parameters")):
        param_type = field_dict["parameter_type"]
        param_id = field_dict["parameter_id"]
        k = field_dict["parameter_name"]
        v = field_dict["parameter_value"]
        params[param_type][param_id][k] = v
    # read libraries
    libraries = collections.OrderedDict()
    for field_dict in read_wksheet(wkbook.sheet_by_name(LIBRARY_ROOT_TAG)):
        # add params to field dict
        myparams = {}
        myparams.update(params["patient"][field_dict["patient_id"]])
        myparams.update(params["sample"][field_dict["sample_id"]])
        myparams.update(params["library"][field_dict["library_id"]])        
        field_dict["params"] = myparams
        # build object
        library = Library(**field_dict)
        # ensure unique ids
        if library.library_id in libraries:
            raise Exception("Found duplicate library id %s" % (library.library_id))
        libraries[library.library_id] = library
    return libraries
