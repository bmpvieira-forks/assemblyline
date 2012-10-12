'''
Created on Aug 7, 2011

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import subprocess

# projects
import assemblyline.rnaseq.config as config
from assemblyline.rnaseq.base import many_up_to_date, detect_format
from assemblyline.rnaseq.libtable import Library, FRAGMENT_LAYOUT_PAIRED

import assemblyline.rnaseq
_pipeline_dir = assemblyline.rnaseq.__path__[0]

def bash_check_retcode(msg="ERROR"):
    return 'if [ $? -gt 0 ]; then echo "%s"; exit $?; fi' % (msg)

def setup_modules_environment(pipeline, server):
    commands = ["source %s" % (server.modules_init_script)]
    for name in pipeline.modules:
        basename = name.split("/")[0]
        commands.append("module rm %s" % (basename))
        commands.append("module add %s" % (name))
    return commands

def concatenate_sequences(input_files, output_file):
    outfmt = detect_format(output_file)
    if outfmt == "txt":
        prefix = output_file
    else:
        prefix = os.path.splitext(output_file)[0]
    commands = []
    for filename in input_files:
        fmt = detect_format(filename)
        if fmt == "gz":
            cmd = 'gzip -cd %s >> %s' % (filename, prefix)
        elif fmt == "bz2":
            cmd = 'bzip2 -cd %s >> %s' % (filename, prefix)
        elif fmt == "txt":
            cmd = 'cat %s >> %s' % (filename, prefix)
        commands.append(cmd)
        commands.append(bash_check_retcode())
    fmt = detect_format(output_file)
    if fmt == "gz":
        commands.append('gzip %s' % (prefix))
    elif fmt == "bz2":
        commands.append('bzip2 %s' % (prefix))
    commands.append(bash_check_retcode())
    return commands

def get_pbs_header(job_name,
                   num_processors=1,
                   node_processors=1,
                   node_memory=4096,
                   pbs_script_lines=None, 
                   working_dir=None, 
                   walltime=None, 
                   pmem=None, 
                   mem=None,
                   deps=None,
                   email=None,
                   stdout_filename=None,
                   stderr_filename=None):
    '''
    job_name: string name of job
    num_processors: number of processors to submit with (cannot be greater than the number of processors per node)
    node_processors: number of cores available per node
    node_memory: amount of memory per node (MB)
    pbs_script_lines: list of PBS directives to be added to the script
    working_dir: the "working directory" of the job (allows scripts to access files using relative pathnames)
    walltime: the walltime passed to qsub
    pmem: amount of memory allocated to each processor of this job in MB
    mem: amount of total memory to split across all processors (overrides pmem setting)
    deps: 'None' if no dependencies, or a python list of job ids
    email: 'None', or string containing codes 'b', 'a', or 'e' describing when to email
    stdout_filename: string filename for storing stdout
    stderr_filename: string filename for storing stderr    
    '''    
    if pbs_script_lines is None:
        pbs_script_lines = []
    if isinstance(deps, basestring):
        deps = [deps]
    # ensure number of processors can fit on node
    num_processors = min(num_processors, node_processors)
    resource_fields = ["nodes=%d:ppn=%d" % (1, num_processors)]    
    # setup memory resources
    if mem is not None:
        if mem > node_memory:
            logging.warning("Job requested more memory than node supports (%dmb > %dmb)" % (mem, node_memory))        
        mem = min(mem, node_memory)
        resource_fields.append("mem=%dmb" % (mem))
    elif pmem is not None:
        max_pmem = float(node_memory) / num_processors
        if pmem > max_pmem:
            logging.warning("Job requested more memory-per-process (pmem) than node supports (%dmb > %dmb)" % (pmem, max_pmem))
        pmem = min(pmem, max_pmem)
        resource_fields.append("pmem=%dmb" % (pmem))
    else:
        pmem = int(round(float(node_memory) / node_processors, 0))
        mem = pmem * num_processors
        resource_fields.append("mem=%dmb" % (mem))
    # set job walltime
    if walltime is not None:
        resource_fields.append("walltime=%s" % (walltime))
    # make PBS script
    lines = ["#!/bin/sh",
             "#PBS -N %s" % job_name,
             "#PBS -l %s" % (",".join(resource_fields)),
             "#PBS -V"]
    if email is not None:
        lines.append("#PBS -m %s" % (email))
    if stdout_filename is None:
        stdout_filename = "/dev/null"
    lines.append("#PBS -o %s" % (stdout_filename))
    if stderr_filename is None:
        stderr_filename = "/dev/null"        
    lines.append("#PBS -e %s" % (stderr_filename))
    if deps is not None:
        lines.append("#PBS -W depend=afterok:%s" % (":".join([d for d in deps])))    
    lines.extend(pbs_script_lines)
    if working_dir is not None: 
        lines.append("cd %s" % (working_dir))
    return lines

def submit_pbs(shell_commands):
    p = subprocess.Popen("qsub", stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    p.stdin.write('\n'.join(shell_commands))
    job_id = p.communicate()[0]
    return job_id.strip()

def submit_nopbs(shell_commands,
                 stdout_filename=None,
                 stderr_filename=None):
    outfh = None
    errfh = None
    if stdout_filename is not None:
        outfh = open(stdout_filename, "w")
    if stderr_filename is not None:
        errfh = open(stderr_filename, "w")
    command_str = ';'.join(shell_commands)
    retcode = subprocess.call(command_str, shell=True, stdout=outfh, stderr=errfh)
    if stdout_filename is not None:
        outfh.close()
    if stderr_filename is not None:
        errfh.close()
    return retcode

def run(library_xml_file, config_xml_file, server_name, num_processors,
        stdout_file, stderr_file, dryrun):
    #
    # read and validate configuration file
    #
    logging.info("Reading pipeline configuration file")
    pipeline = config.PipelineConfig.from_xml(config_xml_file)
    if server_name not in pipeline.servers:
        logging.error("Server %s not found" % (server_name))
    server = pipeline.servers[server_name]
    #
    # read library file and attach to results
    #
    library = Library.from_xml_file(library_xml_file)
    if not library.is_valid():
        logging.error("Library not valid")
        return config.JOB_ERROR
    logging.info("Analyzing library: %s" % (library.library_id)) 
    # check genome
    if library.species not in pipeline.species:
        logging.error("Library %s genome %s not found" % 
                      (library.library_id, library.species))
        return config.JOB_ERROR
    # get genome
    genome = pipeline.species[library.species]
    # setup results    
    results = config.RnaseqResults(library, server.output_dir)
    #
    # build up sequence of commands
    #
    shell_commands = []
    #
    # get pbs header
    #
    if server.pbs:
        pbs_commands = get_pbs_header(job_name=library.library_id,
                                      num_processors=num_processors,
                                      node_processors=server.node_processors,
                                      node_memory=server.node_mem,
                                      pbs_script_lines=server.pbs_script_lines, 
                                      walltime=config.PBS_JOB_WALLTIME, 
                                      mem=config.PBS_JOB_MEM,
                                      stdout_filename=stdout_file,
                                      stderr_filename=stderr_file)
        shell_commands.extend(pbs_commands)
    #
    # setup environment
    #
    msg = "Setting up modules environment"
    logging.info(msg)
    commands = setup_modules_environment(pipeline, server)
    logging.debug("Commands: %s" % (commands))
    shell_commands.extend(commands)
    #
    # validate pipeline configuration
    #
    msg = "Validating pipeline configuration"
    logging.info(msg)
    args = [sys.executable, 
            os.path.join(_pipeline_dir, "validate_pipeline_config.py"),
            config_xml_file, server_name]
    logging.debug("\targs: %s" % (' '.join(map(str, args))))
    command = ' '.join(map(str, args))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode())   
    #
    # create directories
    #
    if not os.path.exists(results.output_dir):
        logging.info("Creating directory: %s" % (results.output_dir))
        shell_commands.append("mkdir -p %s" % (results.output_dir))
        shell_commands.append(bash_check_retcode())   
    if not os.path.exists(results.tmp_dir):
        logging.info("Creating directory: %s" % (results.tmp_dir))
        shell_commands.append("mkdir -p %s" % (results.tmp_dir))
        shell_commands.append(bash_check_retcode())   
    if not os.path.exists(results.log_dir):
        logging.info("Creating directory: %s" % (results.log_dir))
        shell_commands.append("mkdir -p %s" % (results.log_dir))
        shell_commands.append(bash_check_retcode())   
    #
    # copy/concatenate sequences read1
    #
    input_files = library.read1_files
    output_files = [results.copied_fastq_files[0]]
    skip = many_up_to_date(output_files, input_files)
    msg = "Concatenating/copying read1 sequence files"
    if skip:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        commands = concatenate_sequences(library.read1_files, results.copied_fastq_files[0])
        shell_commands.extend(commands)
    #
    # copy/concatenate sequences read2
    #
    if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:    
        input_files = library.read2_files
        output_files = [results.copied_fastq_files[1]]
        skip = many_up_to_date(output_files, input_files)
        msg = "Concatenating/copying read2 sequence files"
        if skip:
            logging.info("[SKIPPED] %s" % msg)
        else:
            logging.info(msg)
            commands = concatenate_sequences(library.read2_files, results.copied_fastq_files[1])
            shell_commands.extend(commands)
    #
    # Run FASTQC quality assessment
    #
    input_files = results.copied_fastq_files
    output_files = results.fastqc_data_files + results.fastqc_report_files
    skip = many_up_to_date(output_files, input_files)
    msg = "Running FASTQC quality assessment"
    if skip:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info("%s" % msg)       
        num_threads = min(num_processors, len(results.copied_fastq_files))
        args = [pipeline.fastqc_bin, "--threads", num_threads, "-o", results.output_dir]
        args.extend(results.copied_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # Filter reads that map to abundant sequences
    #
    input_files = results.fastqc_data_files + results.copied_fastq_files
    output_files = results.filtered_fastq_files + [results.sorted_abundant_bam_file]
    skip = many_up_to_date(output_files, input_files)
    msg = "Filtering reads that map to abundant sequences"
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [sys.executable, 
                os.path.join(_pipeline_dir, "filter_abundant_sequences.py"),
                "--num-processors", num_processors,
                ','.join(results.fastqc_data_files),
                ','.join(results.copied_fastq_files),
                ','.join(results.filtered_fastq_files),
                results.sorted_abundant_bam_file,
                os.path.join(server.references_dir, genome.get_path("abundant_bowtie2_index")),
                pipeline.picard_dir,
                results.tmp_dir]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # determine the fragment size distribution of the reads
    #
    msg = "Estimating fragment size distribution"
    input_files = results.filtered_fastq_files
    output_files = [results.frag_size_dist_file, results.frag_size_dist_plot_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)        
        args = [sys.executable, 
                os.path.join(_pipeline_dir, "estimate_fragment_size_distribution.py"),
                "--bowtie-bin", pipeline.bowtie_bin,
                "--min-fragment-size", pipeline.min_fragment_size,
                "--max-fragment-size", pipeline.max_fragment_size,
                "--default-mean", pipeline.fragment_size_mean_default,
                "--default-stdev", pipeline.fragment_size_stdev_default,
                "-p", num_processors,
                "-n", 1000000,
                os.path.join(server.references_dir, genome.get_path("fragment_size_bowtie_index")),
                results.frag_size_dist_file,
                results.frag_size_dist_plot_file]
        args.extend(results.filtered_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # align reads with tophat
    #
    input_files = results.filtered_fastq_files + [results.frag_size_dist_file]
    output_files = [results.tophat_bam_file]
    msg = "Aligning reads with Tophat"
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info("%s" % (msg))
        args = [sys.executable, os.path.join(_pipeline_dir, "run_tophat.py"),
                "-p", num_processors,
                "--library-type", library.library_type, 
                '--rg-id', library.library_id,
                '--rg-sample="%s"' % library.sample_id,
                '--rg-library="%s"' % library.library_id,
                '--rg-description="%s"' % library.description,
                '--rg-center="%s"' % library.study_id]
        for arg in pipeline.tophat_args:
            # substitute species-specific root directory
            species_arg = arg.replace("${SPECIES}", os.path.join(server.references_dir, genome.root_dir)) 
            args.extend(['--tophat-arg="%s"' % species_arg])
        args.extend([results.tophat_dir,
                     os.path.join(server.references_dir, genome.get_path("genome_bowtie2_index")),
                     results.frag_size_dist_file])
        args.extend(results.filtered_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)        
        shell_commands.append(bash_check_retcode())   
    #
    # run picard diagnostics for alignment results
    #
    msg = "Collecting alignment metrics with Picard"    
    input_files = [results.tophat_bam_file]
    output_files = [results.alignment_summary_metrics,
                    results.quality_by_cycle_metrics,
                    results.quality_distribution_metrics]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "CollectMultipleMetrics.jar"),
                "INPUT=%s" % (results.tophat_bam_file),
                "REFERENCE_SEQUENCE=%s" % os.path.join(server.references_dir, genome.get_path("genome_lexicographical_fasta_file")),
                "OUTPUT=%s" % (os.path.join(results.output_dir, "picard")),
                "ASSUME_SORTED=TRUE",
                "TMP_DIR=%s" % results.tmp_dir,
                "VALIDATION_STRINGENCY=SILENT"]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # run picard for rna-seq diagnostics
    #
    msg = "Collecting RNA-Seq metrics with Picard"
    input_files = [results.tophat_bam_file]
    output_files = [results.rnaseq_metrics]
    skip = many_up_to_date(output_files, input_files)
    if skip:    
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        args = ["java", "-jar", 
                os.path.join(pipeline.picard_dir, "CollectRnaSeqMetrics.jar"),
                "INPUT=%s" % (results.tophat_bam_file),
                "REF_FLAT=%s" % (os.path.join(server.references_dir, genome.get_path("gene_annotation_refflat"))),
                "RIBOSOMAL_INTERVALS=%s" % (os.path.join(server.references_dir, genome.get_path("picard_ribosomal_intervals"))),
                "STRAND_SPECIFICITY=%s" % config.get_picard_strand_specificity(library.library_type),
                "REFERENCE_SEQUENCE=%s" % (os.path.join(server.references_dir, genome.get_path("genome_lexicographical_fasta_file"))),
                "OUTPUT=%s" % (results.rnaseq_metrics),
                "CHART_OUTPUT=%s" % (results.rnaseq_metrics_pdf),
                "TMP_DIR=%s" % results.tmp_dir]
                #"ASSUME_SORTED=TRUE",
                #"VALIDATION_STRINGENCY=SILENT"]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # generate genome coverage bedgraph  
    #
    msg = "Generating coverage bedgraph file"
    input_files = [results.tophat_bam_file]
    output_files = [results.coverage_bedgraph_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [os.path.join(pipeline.bedtools_dir, "genomeCoverageBed"),
                "-bg", "-split", "-g", (os.path.join(server.references_dir, genome.get_path("chrom_sizes"))),
                "-ibam", results.tophat_bam_file, ">", results.coverage_bedgraph_file]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # convert bedgraph to bigwig coverage file
    #
    msg = "Create bigWig file to display coverage"
    input_files = [results.coverage_bedgraph_file]
    output_files = [results.coverage_bigwig_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [os.path.join(pipeline.ucsc_dir, "bedGraphToBigWig"),
                results.coverage_bedgraph_file,
                os.path.join(server.references_dir, genome.get_path("chrom_sizes")),
                results.coverage_bigwig_file]
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # run cufflinks to assembly transcriptome
    #
    msg = "Estimating known transcript abundances with Cufflinks"
    output_files = [results.cufflinks_gtf_file]
    input_files = [results.frag_size_dist_file, results.tophat_bam_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        if not os.path.exists(results.cufflinks_dir):
            logging.info("\tcreating directory: %s" % (results.cufflinks_dir))
            os.makedirs(results.cufflinks_dir)
        args = [sys.executable, os.path.join(_pipeline_dir, "run_cufflinks.py"),
                "--cufflinks-bin", pipeline.cufflinks_bin,
                "-p", num_processors,
                "-L", library.library_id,
                "--library-type", library.library_type]
        if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
            args.append("--learn-frag-size")
        for arg in pipeline.cufflinks_args:
            # substitute species-specific root directory
            species_arg = arg.replace("${SPECIES}", os.path.join(server.references_dir, genome.root_dir)),
            args.append('--cufflinks-arg="%s"' % (species_arg))
        args.extend([results.tophat_bam_file,
                     results.cufflinks_dir,
                     results.frag_size_dist_file])
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # write job finished file
    #
    msg = "Writing job complete file"
    output_files = [results.job_done_file]
    input_files = [results.coverage_bigwig_file,
                   results.tophat_bam_file, 
                   results.cufflinks_gtf_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.info("[SKIPPED] %s" % msg)
    else:
        logging.info(msg)
        command = "touch %s" % (results.job_done_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # cleanup intermediate files
    #
    msg = "Cleaning up tmp files"
    logging.info(msg)
    command = "rm -rf %s" % (results.tmp_dir)
    shell_commands.append(command)    
    shell_commands.append(bash_check_retcode())   
    #
    # show commands
    # 
    for command in shell_commands:
        print command
    #
    # Execute job
    #
    if not dryrun:
        if server.pbs:
            submit_pbs(shell_commands)
        else:
            return submit_nopbs(shell_commands, stdout_file, stderr_file)
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("Assemblyline version %s" % (assemblyline.__version__))
    logging.info("=============================")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("--stdout", dest="stdout_file", default=None)
    parser.add_argument("--stderr", dest="stderr_file", default=None)
    parser.add_argument("--dryrun", dest="dryrun", action="store_true", 
                        default=False)
    parser.add_argument("library_xml_file")
    parser.add_argument("config_xml_file")
    parser.add_argument("server_name")
    args = parser.parse_args()
    return run(args.library_xml_file, args.config_xml_file, 
               args.server_name, args.num_processors,
               args.stdout_file, args.stderr_file,
               args.dryrun)

if __name__ == '__main__': 
    sys.exit(main())

