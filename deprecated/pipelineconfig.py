'''
Created on Feb 16, 2010

@author: mkiyer
'''
import os
import re
import logging

import xml.etree.cElementTree as etree

def read_sample_pools(line_iter):
    # first line is sample names
    sample_names = line_iter.next().strip().split('\t')[1:]
    # next is sample groups
    sample_groups = line_iter.next().strip().split('\t')[1:]
    # next is color
    colors = line_iter.next().strip().split('\t')[1:]
    # next is library_ids
    library_ids = [s.split(',') for s in line_iter.next().strip().split('\t')[1:]]
    # now build property list and dictionary
    propdict = {}
    for row in line_iter:
        fields = row.strip().split('\t')
        prop = fields[0]
        prop_values = fields[1:]
        propdict[prop] = prop_values
    sample_pools = type('SamplePoolInfo', (object,), 
                        dict(names=sample_names,
                             groups=sample_groups,
                             colors=colors,
                             library_ids=library_ids,
                             properties=propdict))
    return sample_pools

def find_file(path, default_path):
    if os.path.exists(path):
        return os.path.abspath(path)        
    hacked_path = fix_path(path)
    if os.path.exists(hacked_path):
        return os.path.abspath(hacked_path)
    # prepend the prefix and repeat search
    joined_path = os.path.join(default_path, path)        
    if os.path.exists(joined_path):
        return os.path.abspath(joined_path)
    hacked_path = fix_path(joined_path)
    if os.path.exists(hacked_path):
        return os.path.abspath(hacked_path)
    return None

def fix_path(path):
    # hacked garbage
    hacked_path = re.sub('archive11', 'data2', path, count=1)
    if os.path.exists(hacked_path):
        return hacked_path
    hacked_path = re.sub('data2', 'archive11', path, count=1)
    if os.path.exists(hacked_path):
        return hacked_path        
    return path

def get_sample_pools(filename):
    logging.debug("Reading sample pools from %s" % filename)
    pools = read_sample_pools(open(filename))
    logging.debug("Number of samples: %d" % len(pools.names))
    logging.debug("Number of libraries: %d" % sum([len(x) for x in pools.library_ids]))
    return pools

def parse_yesno_text(s):
    return True if s.upper() == 'YES' else False


class PipelineConfig(object):
    
    def __init__(self, xml_file):        
        # remember xml file
        self.xml_file = xml_file
        # fill these variables with values from the xml
        self.parse_xml(xml_file)

    def parse_xml(self, xml_file):
        tree = etree.parse(xml_file)
        root = tree.getroot()    
        assert root.tag == 'pipeline'        
        # general options
        output_dir = root.findtext('output_dir')
        self.output_dir = output_dir if output_dir != None else None        
        data_dir = root.findtext('data_dir')
        if (data_dir != None) and (data_dir != ''):
            assert os.path.exists(data_dir)
        self.data_dir = data_dir if data_dir != None else None
        self.genome_version = root.findtext('genome_version')    
            
        # sample tag
        self.sample_pools_file = root.findtext('samples')
        self.sample_pools = get_sample_pools(self.sample_pools_file)

        # point to tuxedo data
        elem = root.find('tuxedo')
        self.tuxedo_root_dir = elem.get('dir')
        self.tuxedo_data_dir = os.path.join(self.tuxedo_root_dir, elem.findtext('data'))

        # cuffcompare tags
        elem = root.find('cuffcompare')
        self.cuffcompare_executable = os.path.expandvars(elem.get('bin'))
        self.cuffcompare_output_path = fix_path(elem.findtext('output_path'))
        elem2 = elem.find('ref_gtf')
        self.cuffcompare_ref_enabled = parse_yesno_text(elem2.get('enabled'))
        if self.cuffcompare_ref_enabled:
            self.cuffcompare_ref_gtf = find_file(elem2.text, data_dir)
        else:
            self.cuffcompare_ref_gtf = None
        self.cuffcompare_analysis_dir = root.findtext('cuffcompare_analysis_dir')
  

if __name__ == '__main__':
    import sys
    p = PipelineConfig()
    p.parse_xml(sys.argv[1])
    print p
