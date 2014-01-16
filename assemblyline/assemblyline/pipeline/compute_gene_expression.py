'''
Created on May 16, 2013

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import xml.etree.cElementTree as etree

from assemblyline.lib.base import Library
import assemblyline.lib
_script_dir = assemblyline.lib.__path__[0]

class Config(object):
    @staticmethod
    def from_xml(xmlfile):
        tree = etree.parse(xmlfile)  
        root = tree.getroot()
        c = Config()
        c.init_script = root.findtext("init_script")
        # pbs
        c.pbs = False
        c.node_mem = None
        c.node_processors = None
        c.pbs_script_lines = []
        pbs_elem = root.find("pbs")
        has_pbs = (pbs_elem.get("use") == "yes")
        if has_pbs:
            c.pbs = True
            c.node_mem = float(pbs_elem.findtext("node_mem"))
            c.node_processors = int(pbs_elem.findtext("node_processors"))
            for line_elem in pbs_elem.findall("script_line"):
                c.pbs_script_lines.append(line_elem.text)
        # modules
        c.modules = []
        modules_elem = root.find("modules")
        for line_elem in modules_elem.findall("module"):
            c.modules.append(line_elem.text)
        return c

def bash_log(msg, level="DEBUG"):
    return 'echo "[`date`] - %s - %s" >&2' % (level, msg)

def bash_check_retcode(msg="ERROR", level="ERROR"):
    echo_command = bash_log(msg, level)
    return 'ERRCODE=$?; if [ $ERRCODE -gt 0 ]; then %s; exit $ERRCODE; fi' % (echo_command)

def bash_check_file_and_exit(filename, msg="File exists. Skipping execution."):
    echo_command = bash_log(msg, level="INFO")
    return 'if [ -e %s ]; then %s; exit 1; fi' % (filename, echo_command)

def get_pbs_header(job_name,
                   node_processors=1,
                   node_memory=4096,
                   pbs_script_lines=None, 
                   working_dir=None, 
                   deps=None,
                   stdout_filename=None,
                   stderr_filename=None):
    '''
    job_name: string name of job
    node_processors: number of cores available per node
    node_memory: amount of memory per node (MB)
    pbs_script_lines: list of PBS directives to be added to the script
    working_dir: the "working directory" of the job (allows scripts to access files using relative pathnames)
    deps: 'None' if no dependencies, or a python list of job ids
    stdout_filename: string filename for storing stdout
    stderr_filename: string filename for storing stderr    
    '''    
    if pbs_script_lines is None:
        pbs_script_lines = []
    if isinstance(deps, basestring):
        deps = [deps]
    # add PBS parameters
    lines = ["#PBS -N %s" % job_name]
    lines.extend(pbs_script_lines)
    if stdout_filename is None:
        stdout_filename = "/dev/null"
    lines.append("#PBS -o %s" % (stdout_filename))
    if stderr_filename is None:
        stderr_filename = "/dev/null"        
    lines.append("#PBS -e %s" % (stderr_filename))
    if deps is not None:
        lines.append("#PBS -W depend=afterok:%s" % (":".join([d for d in deps])))    
    if working_dir is not None: 
        lines.append("cd %s" % (working_dir))
    return lines

def create_shell_commands(lib, config, output_dir,
                          log_prefix,
                          node_mem):
    #
    # build up sequence of commands
    #
    shell_commands = []
    shell_commands.append("#!/bin/sh")
    #
    # get pbs header
    #
    if config.pbs:
        stdout_file = os.path.join(output_dir, "%s.stdout" % (log_prefix))
        stderr_file = os.path.join(output_dir, "%s.stderr" % (log_prefix))
        pbs_commands = get_pbs_header(job_name=lib.library_id,
                                      node_processors=config.node_processors,
                                      node_memory=config.node_mem,
                                      pbs_script_lines=config.pbs_script_lines, 
                                      stdout_filename=stdout_file,
                                      stderr_filename=stderr_file)
        shell_commands.extend(pbs_commands)
    #
    # intro logging messages
    #
    shell_commands.append(bash_log("AssemblyLine RNA-Seq Pipeline", "INFO"))
    shell_commands.append(bash_log("Library '%s'" % (lib.library_id), "INFO"))
    #
    # setup environment
    #
    shell_commands.append("source %s" % (config.init_script))
    shell_commands.append("module purge")
    for name in config.modules:        
        shell_commands.append("module add %s" % (name))
    #
    # create directories
    #
    tmp_dir = os.path.join(output_dir, "tmp")
    msg = "Creating directory: %s" % (tmp_dir)
    shell_commands.append(bash_log(msg, "INFO"))
    shell_commands.append("mkdir -p %s" % (tmp_dir))
    shell_commands.append(bash_check_retcode(msg))
    return shell_commands

def cufflinks(lib, gtf_file, config):
    # setup output files
    output_dir = os.path.join(config.output_dir, lib.library_id)
    tmp_dir = os.path.join(output_dir, "tmp")
    job_done_file = os.path.join(output_dir, "cufflinks.done")
    script_file = os.path.join(output_dir, "cufflinks.sh")
    #COUNTS = FPKM*transcript_length/1000000000 * MASS
    # init shell commands
    shell_commands = create_shell_commands(lib, config, output_dir, 
                                           log_prefix="cufflinks",
                                           node_mem=8000)    
    # check that job is not already done
    shell_commands.append(bash_log("Checking if job already complete"))
    shell_commands.append(bash_check_file_and_exit(job_done_file))
    # run cufflinks for gene expression
    msg = "Counting reads across genes with Cufflinks"        
    shell_commands.append(bash_log(msg, "INFO"))
    log_file = os.path.join(output_dir, 'cufflinks.log')
    args = ["cufflinks",
            "-q",
            "--no-update-check",
            "-p", 1,
            "--output-dir", output_dir,
            "--GTF", gtf_file,
            "--max-frag-multihits", 1,
            "--no-effective-length-correction",
            "--max-bundle-length", 5000000,
            "--max-bundle-frags", 20000000,
            lib.bam_file,
            "> %s 2>&1" % (log_file)]
    command = ' '.join(map(str, args))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))
    # write job done file
    msg = "Writing job done file"
    command = "touch %s" % (job_done_file)    
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))
    # cleanup
    msg = "Cleaning up tmp files"
    if config.keep_tmp:
        shell_commands.append(bash_log("[SKIPPED] %s because '--keep-tmp' flag set" % (msg), "INFO"))
    else:
        shell_commands.append(bash_log(msg, "INFO"))
        command = "rm -rf %s" % (tmp_dir)
        shell_commands.append(command)    
        shell_commands.append(bash_check_retcode())
    # write commands 
    if config.dryrun:
        for command in shell_commands:
            print command
    else:
        f = open(script_file, "w")
        for command in shell_commands:
            print >>f, command
        f.close()

def htseq_count(lib, gtf_file, config):
    # setup output files
    output_dir = os.path.join(config.output_dir, lib.library_id)
    tmp_dir = os.path.join(output_dir, "tmp")
    job_done_file = os.path.join(output_dir, "htseq_count.done")
    script_file = os.path.join(output_dir, "htseq_count.sh")
    # htseq parameters
    run_as_pe = False
    extra_args = ['-m union', '-s no']
    # init shell commands
    shell_commands = create_shell_commands(lib, config, output_dir,
                                           log_prefix="htseq_count") 
    # check that job is not already done
    shell_commands.append(bash_log("Checking if job already complete"))
    shell_commands.append(bash_check_file_and_exit(job_done_file))
    #
    # run htseq-count for gene expression
    #
    msg = "Counting reads across genes with htseq-count"        
    shell_commands.append(bash_log(msg, "INFO"))
    output_file = os.path.join(output_dir, "htseq.txt")
    log_file = os.path.join(output_dir, 'htseq.log')
    args = ["python", os.path.join(_script_dir, "htseq_count.py"),
            "--tmp-dir", tmp_dir]
    if run_as_pe:
        args.append("--pe")
    for arg in extra_args:
        args.append('--arg="%s"' % arg)
    args.extend([gtf_file, lib.bam_file, output_file,
                 '> %s 2>&1' % (log_file)])                     
    command = ' '.join(map(str, args))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))
    # write job done file
    msg = "Writing job done file"
    command = "touch %s" % (job_done_file)    
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))
    #
    # cleanup
    #
    msg = "Cleaning up tmp files"
    if config.keep_tmp:
        shell_commands.append(bash_log("[SKIPPED] %s because '--keep-tmp' flag set" % (msg), "INFO"))
    else:
        shell_commands.append(bash_log(msg, "INFO"))
        command = "rm -rf %s" % (tmp_dir)
        shell_commands.append(command)    
        shell_commands.append(bash_check_retcode())   
    # write commands 
    if config.dryrun:
        for command in shell_commands:
            print command
    else:
        f = open(script_file, "w")
        for command in shell_commands:
            print >>f, command
        f.close()

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("AssemblyLine %s" % (assemblyline.__version__))
    logging.info("----------------------------------")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", dest="mode", choices=['htseq', 'cufflinks'],
                        default='htseq')
    parser.add_argument("--dryrun", dest="dryrun", action="store_true", 
                        default=False)
    parser.add_argument("--keep-tmp", dest="keep_tmp", action="store_true", 
                        default=None)
    parser.add_argument('-o', '--output-dir', dest="output_dir", default=None)
    parser.add_argument("config_xml_file")
    parser.add_argument('library_table')
    parser.add_argument("gtf_file")
    args = parser.parse_args()
    if not os.path.exists(args.config_xml_file):
        parser.error("config xml file '%s' not found" % (args.config_xml_file))
    if not os.path.exists(args.library_table):
        parser.error("library table file '%s' not found" % (args.library_table))
    if not os.path.exists(args.gtf_file):
        parser.error("gtf file '%s' not found" % (args.gtf_file))
    gtf_file = os.path.abspath(args.gtf_file)
    # read configuration
    logging.info("Reading configuration file")
    config = Config.from_xml(args.config_xml_file)
    if args.output_dir is not None:
        config.output_dir = os.path.abspath(args.output_dir)
    else:
        config.output_dir = os.getcwd()
    config.keep_tmp = args.keep_tmp
    config.dryrun = args.dryrun
    # create output directory
    if not os.path.exists(config.output_dir):
        os.makedirs(config.output_dir)
    # read library table
    logging.info("Parsing library table")
    num_libs = 0
    for lib in Library.from_file(args.library_table):
        if not os.path.exists(lib.bam_file):
            logging.warning("\t[SKIPPED] Library %s BAM file not found" % (lib.library_id))
            continue
        # setup library dir
        lib_dir = os.path.join(config.output_dir, lib.library_id)
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        if args.mode == 'htseq':
            htseq_count(lib, gtf_file, config)
        elif args.mode == 'cufflinks':
            cufflinks(lib, gtf_file, config)
        num_libs += 1
    logging.info("Found %d libraries" % (num_libs))
    logging.info("Done")
    return 0

if __name__ == '__main__':
    sys.exit(main())