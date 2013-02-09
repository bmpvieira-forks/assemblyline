'''
Created on Dec 19, 2012

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import urllib

# projects
import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.libtable import Library
from assemblyline.rnaseq.lib.inspect import RnaseqLibraryMetrics
from assemblyline.rnaseq.lib.libtable import FR_UNSTRANDED

import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

STRAND_COLOR_MAP = {".": "0,128,0",
                    "+": "255,0,0",
                    "-": "0,0,255"}

UCSC_BASEURL = "http://genome.ucsc.edu/cgi-bin/hgTracks?"

def get_tinyurl(url):
    apiurl = "http://tinyurl.com/api-create.php?url="
    tinyurl = urllib.urlopen(apiurl + url).read()
    return tinyurl

#http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr21:33038447-33041505&hgct_customText
#=track%20type=bigBed%20name=myBigBedTrack%20description=%22a%20bigBed%20track%22%20visibility=
#full%20bigDataUrl=http://genome.ucsc.edu/goldenPath/help/examples/bigBedExample.bb
#
#def get_custom_track_url(track_string,
#                         pos=None,
#                         ucsc_url=default_ucsc_url,
#                         ucsc_db=default_ucsc_db,
#                         **track_kwargs):
#    # specify UCSC URL parameters
#    url_params = ['db=%s' % (ucsc_db)]
#    if pos is not None:
#        url_params.append('position=%s' % pos)
#    url_params.append('hgct_customText=%s' % (track_string))
#    url_string = ucsc_url + '&'.join(url_params)
#    return url_string
#
#def get_track_string(name, desc, 
#                     external_ip,
#                     track_file,
#                     track_type,
#                     **track_kwargs):
#    # convert to html url format
#    local_url = "%s%s" % (external_ip, urllib.pathname2url(track_file))
#    # specify track line as string
#    track_params = ['track', 
#                    'type=%s' % track_type, 
#                    urllib.quote('name="%s"' % (name)),
#                    urllib.quote('description="%s"' % (desc))]
#    for k,v in track_kwargs.iteritems():
#        track_params.append(urllib.quote('%s=%s' % (k,v)))
#    track_params.append(urllib.quote('bigDataUrl="%s"' % (local_url)))
#    track_string = '%20'.join(track_params)
#    return track_string

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxWindowToDraw", dest="maxWindowToDraw", type=int, default=200000)
    parser.add_argument("--pairSearchRange", dest="pairSearchRange", type=int, default=20000)
    parser.add_argument("--position", dest="position", default=None)
    parser.add_argument("--baseurl", dest="baseurl", default=None)
    parser.add_argument("--scratch-dir", dest="scratch_dir", default=None)
    parser.add_argument("--httplink", dest="httplink", action="store_true", default=False)
    parser.add_argument("--tinyurl", dest="tinyurl", action="store_true", default=False)
    parser.add_argument("library_dir")
    args = parser.parse_args()
    # check command line
    if not os.path.exists(args.library_dir):
        parser.error("library_dir '%s' not found" % (args.library_dir))
    if not os.path.isdir(args.library_dir):
        parser.error("library_dir '%s' not a directory" % (args.library_dir))
    if args.scratch_dir is None:
        scratch_dir = os.getcwd()
    elif not os.path.exists(args.scratch_dir):
        parser.error("scratch_dir '%s' not found" % (args.scratch_dir))
    elif not os.path.isdir(args.scratch_dir):
        parser.error("scratch_dir '%s' not not a directory" % (args.scratch_dir))
    else:
        scratch_dir = os.path.abspath(args.scratch_dir)
    if not args.baseurl:
        parser.error("--baseurl not specified")
    if args.position is not None:
        position = args.position.replace(',', '')
    else:
        position = None
    # read library file
    library_dir = os.path.abspath(args.library_dir)
    library_xml_file = os.path.join(library_dir, config.LIBRARY_XML_FILE)
    library = list(Library.from_xml_file(library_xml_file))[0]
    # read pipeline config file
    pipeline_xml_file = os.path.join(library_dir, config.CONFIG_XML_FILE)
    pipeline = config.PipelineConfig.from_xml(pipeline_xml_file)
    # check genome
    if library.species not in pipeline.genomes:
        logging.error("Library %s genome %s not found" % 
                      (library.library_id, library.species))
        return config.JOB_ERROR
    genome = pipeline.genomes[library.species]
    # get results
    results = config.RnaseqResults(library, pipeline, library_dir)
    # predict library type
    obj = RnaseqLibraryMetrics.from_file(results.library_metrics_file)
    library_type = obj.predict_library_type()
    # build track lines
    track_lines = []
    #
    # Tophat (bowtie2) BAM file
    #
    ucsc_url = "%s%s" % (args.baseurl, results.tophat_bam_file)
    track_name = "tophat2_bam_%s" % (results.library_id)
    track_desc = "Tophat2 BAM for %s" % (results.library_id)
    track_line = ' '.join(['track type=bam',
                           'name="%s"' % (track_name),
                           'description="%s"' % (track_desc),
                           'visibility=hide',
                           'pairEndsByName=.',
                           'pairSearchRange=%d' % (args.pairSearchRange),
                           'bamColorMode=strand',
                           'bamGrayMode=unpaired',
                           'maxWindowToDraw=%s' % (args.maxWindowToDraw),
                           'bigDataUrl="%s"' % (ucsc_url)])
    track_lines.append(track_line)
    #
    # Tophat Fusion (bowtie1) BAM file
    #
    ucsc_url = "%s%s" % (args.baseurl, results.tophat_fusion_bam_file)
    track_name = "tophatfus_bam_%s" % (results.library_id)
    track_desc = "Tophatfus BAM for %s" % (results.library_id)
    track_line = ' '.join(['track type=bam',
                           'name="%s"' % (track_name),
                           'description="%s"' % (track_desc),
                           'visibility=hide',
                           'pairEndsByName=.',
                           'pairSearchRange=%d' % (args.maxWindowToDraw),
                           'bamColorMode=strand',
                           'bamGrayMode=unpaired',
                           'maxWindowToDraw=%s' % (args.maxWindowToDraw),
                           'bigDataUrl="%s"' % (ucsc_url)])
    track_lines.append(track_line)
    #
    # coverage files (bigwig)
    # 
    if library_type == FR_UNSTRANDED:
        strand_files = [(".", results.coverage_bigwig_prefix + ".bw")]
    else:
        strand_files = [("+", results.coverage_bigwig_prefix + "_pos.bw"),
                        ("-", results.coverage_bigwig_prefix + "_neg.bw")] 
    for strand,bigwig_file in strand_files:
        ucsc_url = "%s%s" % (args.baseurl, bigwig_file)
        track_name = "cov_%s_%s" % (results.library_id, strand)
        track_desc = "Coverage for %s strand=[%s] (RPM)" % (results.library_id, strand)
        track_line = ' '.join(['track type=bigWig',
                               'name="%s"' % (track_name),
                               'description="%s"' % (track_desc),
                               'visibility=full',
                               'color=%s' % (STRAND_COLOR_MAP[strand]),
                               'autoScale=on',
                               'maxHeightPixels=64:64:11',
                               'bigDataUrl="%s"' % (ucsc_url)])
        track_lines.append(track_line)
    #
    # junction bigbed file
    #
    ucsc_url = "%s%s" % (args.baseurl, results.junctions_bigbed_file)
    track_name = "junc_%s" % (results.library_id)
    track_desc = "Splice junctions for %s" % (results.library_id)
    track_line = ' '.join(['track type=bigBed',
                           'name="%s"' % (track_name),
                           'description="%s"' % (track_desc),
                           'visibility=pack',
                           'bigDataUrl="%s"' % (ucsc_url)])
    track_lines.append(track_line)
    #
    # variant calls (VCF)
    #
    ucsc_url = "%s%s" % (args.baseurl, results.varscan_snv_bgzip_file)
    track_name = "snv_%s" % (results.library_id)
    track_desc = "SNV for %s" % (results.library_id)
    track_line = ' '.join(['track type=vcfTabix',
                           'name="%s"' % (track_name),
                           'description="%s"' % (track_desc),
                           'visibility=pack',
                           'bigDataUrl="%s"' % (ucsc_url)])
    track_lines.append(track_line)
    #
    # http links
    #
    track_file = os.path.abspath(os.path.join(scratch_dir, "%s.ucsc_tracks.txt" % (results.library_id)))
    track_file_url = "%s%s" % (args.baseurl, urllib.pathname2url(track_file))
    url_params = ['%sdb=%s' % (UCSC_BASEURL, genome.ucsc_db)]
    if position is not None:
        url_params.append('position=%s' % (urllib.quote(position)))
    url_params.append('hgt.customText=%s' % track_file_url)
    track_url = '&'.join(url_params)
    if args.tinyurl:
        track_url = get_tinyurl(track_url)
    #
    # output
    #
    if args.httplink:
        f = open(track_file, 'w')
        for line in track_lines:
            print >>f, line
        f.close()
        print track_url
    else:
        if position is not None:
            print "browser position %s" % (position)
        for line in track_lines:
            print line

if __name__ == '__main__':
    sys.exit(main())