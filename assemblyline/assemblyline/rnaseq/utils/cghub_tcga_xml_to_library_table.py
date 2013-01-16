'''
Created on Dec 2, 2012

@author: mkiyer
'''
import os
import sys
import logging
import argparse
import xml.etree.cElementTree as etree

from assemblyline.rnaseq.lib.libtable import Library, \
    FRAGMENT_LAYOUT_SINGLE, FRAGMENT_LAYOUT_PAIRED, \
    POLYA_TRANSCRIPTOME, FR_UNSTRANDED
from assemblyline.rnaseq.lib.base import indent_xml

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", dest="species", default="human")
    parser.add_argument("--library-type", dest="library_type", default=FR_UNSTRANDED)
    parser.add_argument("--library-protocol", dest="library_protocol", default=POLYA_TRANSCRIPTOME)
    parser.add_argument("--param", dest="param_list", action="append", default=None)
    parser.add_argument("--seq-repo", dest="seq_repo", default="tcga")
    parser.add_argument("--xls", dest="write_xls", action="store_true", default=False)
    parser.add_argument("--xml", dest="write_xml", action="store_true", default=False)
    parser.add_argument("xml_file")
    parser.add_argument("root_dir")
    args = parser.parse_args()
    # parse param list
    default_params = {}
    if args.param_list is not None:
        for param in args.param_list:
            k,v = param.split("=")
            default_params[k] = v
    # read libraries
    libraries = []
    tree = etree.parse(args.xml_file)  
    root = tree.getroot()
    for elem in root.findall("Result"):
        study_id = elem.findtext("study")
        patient_id = elem.findtext("participant_id")
        sample_id = elem.findtext("sample_id")
        aliquot_id = elem.findtext("aliquot_id")
        legacy_sample_id = elem.findtext("legacy_sample_id")
        analysis_id = elem.findtext("analysis_id")
        disease_abbr = elem.findtext("disease_abbr")
        analysis_id = analysis_id
        paired_elem = elem.find("experiment_xml/EXPERIMENT_SET/EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED")
        if paired_elem is not None:
            fragment_layout = FRAGMENT_LAYOUT_PAIRED
        else:
            fragment_layout = FRAGMENT_LAYOUT_SINGLE
        bam_files = []
        # find bam file
        files_elem = elem.find("files")
        for file_elem in files_elem.findall("file"):
            filename = file_elem.findtext("filename")
            if os.path.splitext(filename)[-1] != ".bam":
                logging.error("File %s not a BAM file" % (filename))
                continue
            correct_filesize = int(file_elem.findtext("filesize"))
            subpath = os.path.join(disease_abbr, analysis_id, filename)
            path = os.path.join(args.root_dir, subpath)
            if not os.path.exists(path):
                logging.error("Analysis %s not found" % (analysis_id))
                continue
            filesize = os.path.getsize(path)
            if filesize != correct_filesize:
                logging.error("Analysis %s has incorrect filesize" % (analysis_id))
                continue
            bam_files.append(subpath)
        if len(bam_files) == 0:
            logging.error("Analysis %s has no valid bam files" % (analysis_id))
            continue
        if len(bam_files) > 1:
            logging.error("Analysis %s has multiple valid bam files" % (analysis_id))
            continue
        kwargs = {'study_id': study_id,
                  'cohort_id': disease_abbr,
                  'patient_id': patient_id,
                  'sample_id': sample_id,
                  'library_id': aliquot_id,
                  'description': legacy_sample_id,
                  'species': args.species,
                  'library_protocol': args.library_protocol,
                  'library_type': args.library_type,
                  'fragment_layout': fragment_layout,
                  'seq_repo': args.seq_repo,
                  'read1_files': '',
                  'read2_files': '',
                  'bam_files': bam_files[0]}
        kwargs['params'] = {'cancer_type': disease_abbr,
                            'legacy_sample_id': legacy_sample_id}
        kwargs['params'].update(default_params)
        library = Library(**kwargs)
        libraries.append(library)
    # write
    if args.write_xml:
        root = etree.Element("libraries")
        for library in libraries:
            library.to_xml(root)
        # indent for pretty printing
        indent_xml(root)
        print etree.tostring(root)
    elif args.write_xls:
        pass
    else:
        # print file headers
        header_fields = ["study_id", "cohort_id", "patient_id", "sample_id", 
                         "library_id", "description", "species", 
                         "library_protocol", "library_type", 
                         "fragment_layout", "seq_repo", 
                         "read1_files", "read2_files", "bam_files"]
        print '\t'.join(header_fields)
        print '\t'.join(header_fields)
        fields = [study_id, disease_abbr, patient_id, sample_id, aliquot_id,
                  legacy_sample_id, args.species, args.library_protocol, 
                  args.library_type, fragment_layout, 
                  args.seq_repo, "", "", bam_files[0]]
        print '\t'.join(fields)

if __name__ == '__main__':
    sys.exit(main())
        