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


def htseq_count(libs, gtf_file, config):
    # setup output directories
    output_dir = config.output_dir
    tasks_dir = os.path.join(output_dir, "tasks")
    log_dir = os.path.join(output_dir, "logs")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(tasks_dir):
        os.makedirs(tasks_dir)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    # htseq parameters
    run_as_pe = False
    extra_args = ['-m union', '-s no']
    # parse library table
    for lib in libs:
        #
        # build up sequence of commands
        #
        shell_commands = []
        shell_commands.append("#!/bin/sh")
        #
        # get pbs header
        #
        if config.pbs:
            stdout_file = os.path.join(log_dir, "%s.stdout" % (lib.library_id))
            stderr_file = os.path.join(log_dir, "%s.stderr" % (lib.library_id))
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
        lib_output_dir = os.path.join(output_dir, lib.library_id)
        msg = "Creating directory: %s" % (lib_output_dir)
        shell_commands.append(bash_log(msg, "INFO"))
        shell_commands.append("mkdir -p %s" % (lib_output_dir))
        shell_commands.append(bash_check_retcode(msg))  
        tmp_dir = os.path.join(lib_output_dir, "tmp")
        msg = "Creating directory: %s" % (tmp_dir)
        shell_commands.append(bash_log(msg, "INFO"))
        shell_commands.append("mkdir -p %s" % (tmp_dir))
        shell_commands.append(bash_check_retcode(msg))
        #
        # run htseq-count for gene expression
        #
        msg = "Counting reads across genes with htseq-count"        
        shell_commands.append(bash_log(msg, "INFO"))
        output_file = os.path.join(lib_output_dir, "htseq.txt")
        log_file = os.path.join(lib_output_dir, 'htseq.log')
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
        #
        # write commands 
        #
        if config.dryrun:
            for command in shell_commands:
                print command
        else:
            if config.pbs:
                ext = ".pbs"
            else:
                ext = ".sh"
            filename = os.path.join(tasks_dir, "%s%s" % (lib.library_id, ext))
            f = open(filename, "w")
            for command in shell_commands:
                print >>f, command
            f.close()
    return 0

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
    # read configuration
    logging.info("Reading configuration file")
    config = Config.from_xml(args.config_xml_file)
    if args.output_dir is not None:
        config.output_dir = os.path.abspath(args.output_dir)
    else:
        config.output_dir = os.getcwd()
    config.keep_tmp = args.keep_tmp
    config.dryrun = args.dryrun
    # read library table
    logging.info("Parsing library table")
    libs = []
    for lib in Library.from_file(args.library_table):
        if not os.path.exists(lib.bam_file):
            logging.warning("\t[SKIPPED] Library %s BAM file not found" % (lib.library_id))
            continue
        libs.append(lib)
    logging.info("\tfound %d libraries" % (len(libs)))
    if args.mode == 'htseq':
        retcode = htseq_count(libs, args.gtf_file, config)
    elif args.mode == 'cufflinks':
        pass
    logging.info("Done")
    return retcode

if __name__ == '__main__':
    sys.exit(main())