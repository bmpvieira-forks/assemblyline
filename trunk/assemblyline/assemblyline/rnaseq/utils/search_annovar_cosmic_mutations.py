'''
Created on Feb 5, 2013

@author: mkiyer
'''
import os
import sys
import logging
import argparse
import glob

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("mutation_file")
    parser.add_argument("root_dir")
    args = parser.parse_args()
    # read mutations    
    fileh = open(args.mutation_file)
    header = fileh.next()
    mutdict = {}
    for line in fileh:
        fields = line.strip().split('\t')
        mutdict[fields[0]] = fields
    fileh.close()
    # read annovar files
    for library_dir in glob.glob(os.path.join(args.root_dir, "*")):
        if not os.path.isdir(library_dir):
            continue
        library_id = os.path.basename(library_dir)
        annovar_cosmic_file = os.path.join(library_dir, "annovar.hg19_cosmic61_dropped")
        fileh = open(annovar_cosmic_file)
        for line in fileh:
            fields = line.strip().split('\t')
            param_fields = fields[1].split(';')
            params = {}
            for p in param_fields:
                k,valstr = p.split('=')
                vals = valstr.split(',')
                params[k] = vals
            for cosmic_id in params['ID']:
                if cosmic_id in mutdict:
                    output_fields = [library_id] + fields + mutdict[cosmic_id]
                    print '\t'.join(output_fields)
        fileh.close()
    

if __name__ == '__main__':
    sys.exit(main())
