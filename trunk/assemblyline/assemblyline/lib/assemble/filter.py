'''
Created on Aug 21, 2012

@author: mkiyer
'''
import logging
from assemblyline.lib.transcript import NO_STRAND

def filter_transcripts(transcripts, min_length=200):
    """
    removes invalid transcripts
    
    returns new list of transcripts
    """
    # filter transcripts
    new_transcripts = []
    unstranded_multiple_exons = 0
    zero_density = 0
    no_exons = 0
    too_short = 0
    for t in transcripts:
        # check that transcript has exons
        if len(t.exons) == 0:
            no_exons += 1
            continue 
        # check that unstranded transcripts have only one exon
        if t.strand == NO_STRAND:
            if len(t.exons) > 1:
                #logging.debug("Skipping unstranded transcript with multiple exons")
                unstranded_multiple_exons += 1
                continue
        # check length requirement
        if t.length < min_length:
            too_short += 1
            continue
        # check that transcript has positive density
        if t.density > 0:
            new_transcripts.append(t)
        else:
            #logging.debug("Skipping transcript %s with density %f" % (t.id, t.density))
            zero_density += 1
    if no_exons > 0:
        logging.debug("\t\tSkipped %d transcripts with no exons" % (no_exons))
    if unstranded_multiple_exons > 0:
        logging.debug("\t\tSkipped %d unstranded transcripts with multiple exons" % (unstranded_multiple_exons))
    if too_short > 0:
        logging.debug("\t\tSkipped %d transcripts with length < %d" % (too_short, min_length))
    if zero_density > 0:
        logging.debug("\t\tSkipped %d transcripts with zero density" % (zero_density))
    return new_transcripts