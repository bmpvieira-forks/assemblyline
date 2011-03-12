'''
Created on May 4, 2010

@author: mkiyer
'''
import logging
import sys
import xml.etree.cElementTree as etree

# packages

# project imports

# local imports
import cuffcompare

def parse_config_file(xml_file):
    params = {}    
    tree = etree.parse(xml_file)
    root = tree.getroot()
    # data to collect
    elem = root.find('filter')
    params['exempt_exons'] = int(elem.findtext('exempt_exons'))
    logging.info("Transcripts with >= %d exons are exempt" % params['exempt_exons'])
    #params['recurrence'] = int(elem.findtext('recurrence'))
    #logging.info("Exons must appear in >= %d samples" % params['recurrence'])
    params['unique_longest'] = int(elem.findtext('unique_longest'))
    logging.info("A unique genomic region >= %d bp must be present" % params['unique_longest'])
    return params

def filter_transcripts(infileh, outfileh, config_file):
    # read relevant parts in config file
    params = parse_config_file(config_file)
    exempt_exons = params['exempt_exons']
    min_unique_longest = params['unique_longest']
    #trees = build_annotation_interval_trees(params['annotations'])
    #print '\t'.join(["annotated", "length", "num_exons", "recurrence", "unique", "fpkm_lo", "fpkm", "fpkm_hi"])
    passed_features = 0
    failed_features = 0
    passed_ann = [0, 0]
    failed_ann = [0, 0]
    for transcript_features in \
        cuffcompare.parse_gtf_transcript_features(infileh, 
                                                  cuffcompare.classify_gtf_attrs):
        # use a combination of heuristic methods and 
        # results from classification algorithms
        if len(transcript_features) >= exempt_exons:
            new_features = transcript_features
        else:
            new_features = []
            exon_number = 1
            prediction = (max(f.attrs['pred'] for f in transcript_features) > 0)            
            for f in transcript_features:
                unique_longest = f.attrs['ulong']
                if (prediction == True) and (unique_longest >= min_unique_longest):
                    f.attrs['exon_number'] = exon_number
                    exon_number += 1
                    new_features.append(f)
                else:
                    logging.debug("Failed: %s" % (f))
                    failed_features += 1
                    failed_ann[f.attrs['ann']] += 1
        for f in new_features:
            passed_ann[f.attrs['ann']] += 1
            print >>outfileh, f
        passed_features += len(new_features)
    logging.debug("Passed: %d" % passed_features)
    logging.debug("Passed Unannotated: %d" % passed_ann[0])
    logging.debug("Passed Annotated: %d" % passed_ann[1])
    logging.debug("Failed: %d" % failed_features)
    logging.debug("Failed Unannotated: %d" % failed_ann[0])
    logging.debug("Failed Annotated: %d" % failed_ann[1])

if __name__ == '__main__':
    from optparse import OptionParser
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        level=logging.DEBUG)
    optionparser = OptionParser("usage: %prog [options] <config.xml>")
    (options, args) = optionparser.parse_args()
    input_file = args[0]
    config_file = args[1]
    outfileh = sys.stdout    
    filter_transcripts(open(input_file), outfileh, config_file)