'''
Created on Aug 21, 2012

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
import logging
from assemblyline.lib.transcript import NO_STRAND
from assemblyline.lib.base import GTFAttr, FLOAT_PRECISION

def filter_transcripts(transcripts, min_length=250, guided=False):
    """
    removes invalid transcripts
    
    returns new list of transcripts
    """
    # filter transcripts
    new_transcripts = []
    unstranded_multiple_exons = 0
    zero_score = 0
    no_exons = 0
    too_short = 0
    ref_not_guided = 0
    for t in transcripts:
        # determine if this is a reference transcript
        is_ref = bool(int(t.attrs.get(GTFAttr.REF, "0")))
        if is_ref and (not guided):
            ref_not_guided += 1
            continue
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
        # check that transcript has positive score
        if t.score < FLOAT_PRECISION:
            #logging.debug("Skipping transcript %s with score %f" % (t.id, t.score))
            zero_score += 1
            continue
        new_transcripts.append(t)
    if no_exons > 0:
        logging.debug("\t\tSkipped %d transcripts with no exons" % (no_exons))
    if unstranded_multiple_exons > 0:
        logging.debug("\t\tSkipped %d unstranded transcripts with multiple exons" % (unstranded_multiple_exons))
    if too_short > 0:
        logging.debug("\t\tSkipped %d transcripts with length < %d" % (too_short, min_length))
    if zero_score > 0:
        logging.debug("\t\tSkipped %d transcripts with zero score" % (zero_score))
    if ref_not_guided > 0:
        logging.debug("\t\tSkipped %d reference transcripts (not in guided mode)" % (ref_not_guided))
    return new_transcripts