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
import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.base import many_up_to_date, detect_format
from assemblyline.rnaseq.lib.libtable import Library, FRAGMENT_LAYOUT_PAIRED, FRAGMENT_LAYOUT_SINGLE

import assemblyline.rnaseq.pipeline
_pipeline_dir = assemblyline.rnaseq.pipeline.__path__[0]

def bash_log(msg, level="DEBUG"):
    return 'echo "[`date`] - %s - %s" >&2' % (level, msg)

def bash_check_retcode(msg="ERROR", level="ERROR"):
    echo_command = bash_log(msg, level)
    return 'ERRCODE=$?; if [ $ERRCODE -gt 0 ]; then %s; exit $ERRCODE; fi' % (echo_command)

def bash_remove_files(filenames):
    commands = []
    for f in filenames:        
        command = "test -e %s && rm -rf %s" % (f,f)
        commands.append(command)
        commands.append(bash_check_retcode())
    return commands

def bash_set_tmpdir(tmp_dir):
    return 'export TMPDIR=%s' % (tmp_dir)

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
    # add PBS parameters
    lines = ["#PBS -N %s" % job_name,
             "#PBS -l %s" % (",".join(resource_fields))]
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
    if ((stdout_filename is not None) and
        (stdout_filename == stderr_filename)):
        outfh = open(stdout_filename, "w")
        errfh = outfh
    else:
        if stdout_filename is not None:
            outfh = open(stdout_filename, "w")
        if stderr_filename is not None:
            errfh = open(stderr_filename, "w")
    command_str = ';'.join(shell_commands)
    retcode = subprocess.call(command_str, shell=True, stdout=outfh, stderr=errfh)
    if (outfh is not None) and (outfh == errfh):
        outfh.close()
    else:
        if stdout_filename is not None:
            outfh.close()
        if stderr_filename is not None:
            errfh.close()
    return retcode

def run(library_xml_file, config_xml_file, server_name, num_processors,
        stdout_file, stderr_file, dryrun, keep_tmp):
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
    if library.species not in pipeline.genomes:
        logging.error("Library %s genome %s not found" % 
                      (library.library_id, library.species))
        return config.JOB_ERROR
    # resolve genome paths now that server is known
    genome_local = pipeline.genomes[library.species]
    genome_static = genome_local.resolve_paths(server.references_dir)    
    # setup results    
    output_dir = os.path.join(server.output_dir, library.library_id)
    results = config.RnaseqResults(library, output_dir)
    #
    # build up sequence of commands
    #
    shell_commands = []
    #
    # get pbs header
    #
    if server.pbs:
        shell_commands.append("#!/bin/sh")
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
    msg = "Setting TMPDIR environment variable to %s" % (results.tmp_dir)
    logging.debug(msg)
    shell_commands.append(bash_set_tmpdir(results.tmp_dir))
    msg = "Setting up modules environment"
    logging.debug(msg)
    commands = setup_modules_environment(pipeline, server)
    logging.debug("Commands: %s" % (commands))
    shell_commands.extend(commands)
    #
    # validate pipeline configuration
    #
    msg = "Validating pipeline configuration"
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    args = ["python", 
            os.path.join(_pipeline_dir, "validate_pipeline_config.py"),
            config_xml_file, server_name]
    command = ' '.join(map(str, args))
    logging.debug("\tcommand: %s" % (command))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))   
    #
    # create directories
    #
    if not os.path.exists(results.output_dir):
        logging.debug("Creating directory: %s" % (results.output_dir))
        shell_commands.append("mkdir -p %s" % (results.output_dir))
        shell_commands.append(bash_check_retcode())   
    if not os.path.exists(results.tmp_dir):
        logging.debug("Creating directory: %s" % (results.tmp_dir))
        shell_commands.append("mkdir -p %s" % (results.tmp_dir))
        shell_commands.append(bash_check_retcode())   
    if not os.path.exists(results.log_dir):
        logging.debug("Creating directory: %s" % (results.log_dir))
        shell_commands.append("mkdir -p %s" % (results.log_dir))
        shell_commands.append(bash_check_retcode())
    #
    # copy xml files
    #
    msg = "Copying configuration files"    
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    shell_commands.append("cp %s %s" % (library_xml_file, results.library_xml_file))
    shell_commands.append(bash_check_retcode())
    shell_commands.append("cp %s %s" % (config_xml_file, results.config_xml_file))
    shell_commands.append(bash_check_retcode())
    #
    # extract fastq from bam input files
    #
    input_files = library.bam_files
    output_files = results.bam_read1_files + results.bam_read2_files
    skip = many_up_to_date(output_files, input_files)
    msg = "Converting BAM input files to FASTQ"
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        # convert bam to fastq
        for i,prefix in enumerate(results.bam_fastq_prefixes):
            args = ["python",
                    os.path.join(_pipeline_dir, "extract_fastq_from_bam.py"),
                    "--readnum-in-qname",
                    "--tmp-dir", results.tmp_dir,
                    "--fragment-layout", library.fragment_layout,
                    library.bam_files[i], prefix] 
            command = ' '.join(map(str, args))
            shell_commands.append(command)
            shell_commands.append(bash_check_retcode())
    #
    # copy/concatenate sequences read1
    #
    input_files = library.read1_files + results.bam_read1_files
    output_files = [results.copied_fastq_files[0]]
    skip = many_up_to_date(output_files, input_files)
    msg = "Concatenating/copying read1 sequence files"
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        commands = concatenate_sequences(library.read1_files + results.bam_read1_files, 
                                         results.copied_fastq_files[0])
        shell_commands.extend(commands)
    #
    # copy/concatenate sequences read2
    #
    if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:    
        input_files = library.read2_files + results.bam_read2_files
        output_files = [results.copied_fastq_files[1]]
        skip = many_up_to_date(output_files, input_files)
        msg = "Concatenating/copying read2 sequence files"
        if skip:
            logging.debug("[SKIPPED] %s" % msg)
            shell_commands.append(bash_log(msg, "SKIPPED"))
        else:
            logging.debug(msg)
            shell_commands.append(bash_log(msg, "INFO"))
            commands = concatenate_sequences(library.read2_files + results.bam_read2_files, 
                                             results.copied_fastq_files[1])
            shell_commands.extend(commands)
    #
    # Run FASTQC quality assessment
    #
    input_files = results.copied_fastq_files
    output_files = results.fastqc_data_files + results.fastqc_report_files
    skip = many_up_to_date(output_files, input_files)
    msg = "Running FASTQC quality assessment"
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug("%s" % msg)    
        shell_commands.append(bash_log(msg, "INFO"))
        num_threads = min(num_processors, len(results.copied_fastq_files))
        args = ['fastqc', "--threads", num_threads, "-o", results.output_dir]
        args.extend(results.copied_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'fastqc.log')
        command += ' > %s 2>&1' % (log_file)
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
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", 
                os.path.join(_pipeline_dir, "filter_abundant_sequences.py"),
                "--num-processors", num_processors,
                ','.join(results.fastqc_data_files),
                ','.join(results.copied_fastq_files),
                ','.join(results.filtered_fastq_files),
                results.sorted_abundant_bam_file,
                genome_static.abundant_bowtie2_index,
                results.tmp_dir]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'filter_abundant_sequences.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # Count abundant reads
    #
    input_files = [results.sorted_abundant_bam_file]
    output_files = [results.abundant_counts_file]
    skip = many_up_to_date(output_files, input_files)
    msg = "Counting abundant reads"
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'abundant_counts.log')
        command = ("samtools idxstats %s > %s 2> %s" % 
                   (results.sorted_abundant_bam_file,
                    results.abundant_counts_file,
                    log_file))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())    
    #
    # inspect the library and determine characteristics
    #
    msg = "Inspecting library"
    input_files = results.filtered_fastq_files
    output_files = [results.library_metrics_file, 
                    results.frag_size_dist_plot_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg) 
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python",
                os.path.join(_pipeline_dir, "inspect_library.py"),
                "--min-frag-size", pipeline.min_fragment_size,
                "--max-frag-size", pipeline.max_fragment_size,
                "--frag-size-mean", pipeline.fragment_size_mean_default,
                "--frag-size-stdev", pipeline.fragment_size_stdev_default,
                "-n", pipeline.max_inspect_samples,
                "-p", num_processors,
                genome_static.fragment_size_bowtie1_index,
                results.library_metrics_file,
                results.frag_size_dist_plot_file]
        args.extend(results.filtered_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'inspect_library.log')
        command += ' > %s 2>&1' % (log_file)        
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # align reads with tophat fusion
    #
    input_files = results.filtered_fastq_files + [results.library_metrics_file]
    output_files = [results.tophat_fusion_bam_file]
    msg = "Aligning reads with Tophat-Fusion"
    skip = ((not pipeline.tophat_fusion_run) or 
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug("%s" % (msg))
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_tophat.py"),
                "-p", num_processors,
                "--library-type", library.library_type, 
                '--rg-id', library.library_id,
                '--rg-sample="%s"' % library.sample_id,
                '--rg-library="%s"' % library.library_id,
                '--rg-description="%s"' % library.description,
                '--rg-center="%s"' % library.study_id]
        # resolve genome-specific args
        for arg in pipeline.tophat_fusion_args:
            args.append('--tophat-arg="%s"' % genome_static.resolve_arg(arg))
        args.extend([results.tophat_fusion_dir,
                     genome_static.genome_bowtie1_index,
                     results.library_metrics_file])
        args.extend(results.filtered_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'tophat_fusion.log')
        command += ' > %s 2>&1' % (log_file)        
        shell_commands.append(command)        
        shell_commands.append(bash_check_retcode())
    #
    # index tophat fusion bam file
    #
    input_files = [results.tophat_fusion_bam_file]
    output_files = [results.tophat_fusion_bam_index_file]
    msg = "Indexing Tophat-Fusion BAM file"
    skip = ((not pipeline.tophat_fusion_run) or
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug("%s" % (msg))
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["samtools", "index", results.tophat_fusion_bam_file]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)        
        shell_commands.append(bash_check_retcode())        
    #
    # tophat fusion post processing script
    #
    input_files = [results.tophat_fusion_file]
    output_files = [results.tophat_fusion_post_result_file]
    msg = "Post-processing Tophat-Fusion results"
    skip = (((not pipeline.tophat_fusion_run) and
             (not pipeline.tophat_fusion_post_run)) or 
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug("%s" % (msg))
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "tophat_fusion_post.py"),
                "-p", num_processors,
                "-o", results.tophat_fusion_dir]
        args.extend(pipeline.tophat_fusion_post_args)
        args.extend([genome_static.genome_bowtie1_index,
                     genome_static.gene_annotation_refgene,
                     genome_static.gene_annotation_ensgene,
                     results.tophat_fusion_dir,
                     results.library_id])
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'tophat_fusion_post.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)        
        shell_commands.append(bash_check_retcode())
        # cleanup tophat fusion garbage tmp files
        shell_commands.extend(bash_remove_files(results.tophat_fusion_tmp_files))
    #
    # align reads with tophat
    #
    input_files = results.filtered_fastq_files + [results.library_metrics_file]
    output_files = [results.tophat_bam_file]
    msg = "Aligning reads with Tophat"
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug("%s" % (msg))
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_tophat.py"),
                "-p", num_processors,
                "--library-type", library.library_type, 
                '--rg-id', library.library_id,
                '--rg-sample="%s"' % library.sample_id,
                '--rg-library="%s"' % library.library_id,
                '--rg-description="%s"' % library.description,
                '--rg-center="%s"' % library.study_id]
        # resolve genome-specific args
        for arg in pipeline.tophat_args:
            args.append('--tophat-arg="%s"' % genome_static.resolve_arg(arg))
        args.extend([results.tophat_dir,
                     genome_static.genome_bowtie2_index,
                     results.library_metrics_file])
        args.extend(results.filtered_fastq_files)
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'tophat.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)        
        shell_commands.append(bash_check_retcode())
    #
    # index tophat bam file
    #
    input_files = [results.tophat_bam_file]
    output_files = [results.tophat_bam_index_file]
    msg = "Indexing Tophat BAM file"
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug("%s" % (msg))
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["samtools", "index", results.tophat_bam_file]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)        
        shell_commands.append(bash_check_retcode(msg))
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
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["java", "-Xmx4g", "-jar", 
                "$PICARDPATH/CollectMultipleMetrics.jar",
                "INPUT=%s" % (results.tophat_bam_file),
                "REFERENCE_SEQUENCE=%s" % genome_static.genome_lexicographical_fasta_file,
                "OUTPUT=%s" % (os.path.join(results.output_dir, "picard")),
                "ASSUME_SORTED=TRUE",
                "TMP_DIR=%s" % results.tmp_dir,
                "VALIDATION_STRINGENCY=SILENT"]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'picard_collect_multiple_metrics.log')
        command += ' > %s 2>&1' % (log_file)        
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
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python",
                os.path.join(_pipeline_dir, "picard_metrics.py"),
                "--tmp-dir", results.tmp_dir,
                "--picard-dir", "$PICARDPATH",
                results.tophat_bam_file,
                results.library_metrics_file,
                genome_static.gene_annotation_refflat,
                genome_static.picard_ribosomal_intervals,
                genome_static.genome_lexicographical_fasta_file,   
                results.rnaseq_metrics,
                results.rnaseq_metrics_pdf]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'picard_rnaseq_metrics.log')
        command += ' > %s 2>&1' % (log_file)                
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # extract unmapped paired fastq files for realignment
    #
    msg = "Converting unmapped BAM to paired FASTQ for realignment"
    input_files = [results.tophat_unmapped_bam_file]
    output_files = results.unmapped_paired_fastq_files
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        # convert bam to paired fastq
        args = ["python",
                os.path.join(_pipeline_dir, "extract_fastq_from_bam.py"),
                "--fragment-layout", FRAGMENT_LAYOUT_PAIRED,
                "--tmp-dir", results.tmp_dir,
                results.tophat_unmapped_bam_file,
                results.unmapped_paired_fastq_prefix]
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode(msg))
    #
    # map pairs to pathogen sequences
    #
    msg = "Aligning pathogen sequences"
    input_files = results.unmapped_paired_fastq_files
    output_files = [results.pathogen_bam_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg) 
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'bowtie2_align_pathogens.log')        
        args = ["python",
                os.path.join(_pipeline_dir, "bowtie2_paired_align.py"),
                "-p", num_processors,
                "--tmp-dir", results.tmp_dir,
                "--max-frag-size", pipeline.max_fragment_size]
        # extra args
        for arg in pipeline.pathogen_screen_bt2_args:
            args.append('--extra-arg="%s"' % arg)
        args.extend([genome_static.pathogen_bowtie2_index,                     
                     results.pathogen_bam_file])
        args.extend(results.unmapped_paired_fastq_files)
        args.append('> %s 2>&1' % (log_file))
        command = ' '.join(map(str, args))
        logging.debug("\tcommand: %s" % (command))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # count pathogen sequences
    #
    msg = "Counting pathogen sequences"
    input_files = [results.pathogen_bam_file]
    output_files = [results.pathogen_counts_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'pathogen_counts.log')
        command = ("samtools idxstats %s > %s 2> %s" % 
                   (results.pathogen_bam_file,
                    results.pathogen_counts_file,
                    log_file))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())            
    #
    # extract unmapped unpaired fastq files for realignment
    #
    msg = "Converting unmapped BAM to unpaired FASTQ for realignment"
    input_files = [results.tophat_unmapped_bam_file]
    output_files = results.unmapped_unpaired_fastq_files
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        # convert bam to paired fastq
        args = ["python",
                os.path.join(_pipeline_dir, "extract_fastq_from_bam.py"),
                "--fragment-layout", FRAGMENT_LAYOUT_SINGLE,
                "--tmp-dir", results.tmp_dir,
                results.tophat_unmapped_bam_file,
                results.unmapped_unpaired_fastq_prefix]
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode(msg))
    #
    # analyze repeat element content
    #
    msg = "Aligning and counting repeat element sequences"
    input_files = results.unmapped_unpaired_fastq_files
    output_files = [results.repeat_element_counts_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg) 
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'repeat_elements.log')
        args = ["python",
                os.path.join(_pipeline_dir, "bowtie2_unpaired_align.py"),
                "-p", num_processors,
                "--tmp-dir", results.tmp_dir]
        # extra args
        for arg in pipeline.repeat_elements_bt2_args:
            args.append('--extra-arg="%s"' % arg)
        args.extend([genome_static.repbase_bowtie2_index,
                     results.repeat_element_counts_file])
        args.extend(results.unmapped_unpaired_fastq_files)
        args.append('> %s 2>&1' % (log_file))
        command = ' '.join(map(str, args))
        logging.debug("\tcommand: %s" % (command))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode(msg))
    #
    # generate genome coverage maps  
    #
    msg = "Generating coverage maps"
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    args = ["python",
            os.path.join(_pipeline_dir, "coverage_maps.py"),
            "--scale",
            "--tmp-dir", results.tmp_dir,
            results.library_metrics_file,
            results.alignment_summary_metrics,
            genome_static.chrom_sizes,
            results.tophat_bam_file,
            results.coverage_bigwig_prefix]
    logging.debug("\targs: %s" % (' '.join(map(str, args))))
    command = ' '.join(map(str, args))
    log_file = os.path.join(results.log_dir, 'coverage_maps.log')
    command += ' > %s 2>&1' % (log_file)                
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))        
    #
    # generate splice junction bigbed file  
    #
    msg = "Generating splice junction bigbed file"
    input_files = [results.tophat_juncs_file] 
    output_files = [results.junctions_bigbed_file]    
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'bed_to_bigbed.log')
        args = ["python",
                os.path.join(_pipeline_dir, "bed_to_bigbed.py"),
                "--score-to-name",
                "--tmp-dir", results.tmp_dir,
                results.tophat_juncs_file,
                genome_static.chrom_sizes,
                results.junctions_bigbed_file,
                '> %s 2>&1' % (log_file)]
        command = ' '.join(map(str, args))
        logging.debug("\tcommand: %s" % (command))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode()) 
    #
    # run cufflinks to assemble transcriptome
    #
    msg = "Assembling transcriptome with Cufflinks"
    input_files = [results.library_metrics_file, 
                   results.tophat_bam_file]
    output_files = [results.cufflinks_ab_initio_gtf_file]    
    skip = ((not pipeline.cufflinks_ab_initio_run) or
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_cufflinks.py"),
                "-p", num_processors,
                "-L", library.library_id,
                "--library-type", library.library_type]
        if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
            args.append("--learn-frag-size")
        # resolve genome-specific args
        for arg in pipeline.cufflinks_ab_initio_args:
            args.append('--cufflinks-arg="%s"' % genome_static.resolve_arg(arg))
        args.extend([results.tophat_bam_file,
                     results.cufflinks_ab_initio_dir,
                     results.library_metrics_file])
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'cufflinks_ab_initio.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # cufflinks on known genes
    #
    msg = "Estimating known gene expression with cufflinks"
    input_files = [results.library_metrics_file, 
                   results.tophat_bam_file]
    output_files = [results.cufflinks_known_gtf_file]    
    skip = ((not pipeline.cufflinks_known_run) or
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_cufflinks.py"),
                "-p", num_processors,
                "-L", library.library_id,
                "--library-type", library.library_type]
        if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
            args.append("--learn-frag-size")
        # resolve genome-specific args
        for arg in pipeline.cufflinks_known_args:
            args.append('--cufflinks-arg="%s"' % genome_static.resolve_arg(arg))
        args.extend([results.tophat_bam_file,
                     results.cufflinks_known_dir,
                     results.library_metrics_file])
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'cufflinks_known.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # run htseq-count for gene expression
    #
    msg = "Counting reads across genes with htseq-count"
    output_files = [results.htseq_count_known_file]
    input_files = [results.tophat_bam_file]
    skip = ((not pipeline.htseq_count_run) or
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_htseq_count.py"),
                "--tmp-dir", results.tmp_dir]
        if pipeline.htseq_count_pe:
            args.append("--pe")
        for arg in pipeline.htseq_count_args:
            args.append('--arg="%s"' % genome_static.resolve_arg(arg))
        args.extend([genome_static.known_genes_gtf,
                     results.tophat_bam_file,
                     results.htseq_count_known_file])
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'htseq_count_known.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # run picard to remove duplicates
    #
    msg = "Removing duplicate reads with Picard"
    input_files = [results.tophat_bam_file]
    output_files = [results.tophat_rmdup_bam_file]
    skip = ((not pipeline.varscan_run) or
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["java", "-Xmx4g", "-jar", 
                "$PICARDPATH/MarkDuplicates.jar",
                "INPUT=%s" % (results.tophat_bam_file),
                "OUTPUT=%s" % (results.tophat_rmdup_bam_file),
                "METRICS_FILE=%s" % (results.duplicate_metrics),
                "REMOVE_DUPLICATES=true",
                "TMP_DIR=%s" % results.tmp_dir]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'picard_mark_duplicates.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())   
    #
    # run varscan for variant calling
    #
    msg = "Calling variants with VarScan"
    output_files = [results.varscan_snv_file]
    input_files = [results.tophat_rmdup_bam_file]
    skip = ((not pipeline.varscan_run) or
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_varscan.py")]
        for arg in pipeline.varscan_args:
            args.append('--varscan-arg="%s"' % (arg))
        args.extend(["$VARSCANPATH/VarScan.jar", "mpileup2snp",
                     genome_static.genome_lexicographical_fasta_file,
                     results.tophat_rmdup_bam_file,
                     results.varscan_snv_file])
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        log_file = os.path.join(results.log_dir, 'varscan.log')
        command += ' > %s 2>&1' % (log_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # index varscan vcf file
    #
    msg = "Indexing VCF files"
    output_files = [results.varscan_snv_bgzip_file,
                    results.varscan_snv_tabix_file]
    input_files = [results.varscan_snv_file]
    skip = ((not pipeline.varscan_run) or
            many_up_to_date(output_files, input_files))
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        command = 'bgzip -c %s > %s' % (results.varscan_snv_file, 
                                        results.varscan_snv_bgzip_file)
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
        command = 'tabix -p vcf %s' % (results.varscan_snv_bgzip_file) 
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # write job finished file
    #
    msg = "Validating results"
    output_files = [results.job_done_file]
    input_files = [results.library_xml_file, results.config_xml_file]
    skip = many_up_to_date(output_files, input_files)
    if skip:
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", 
                os.path.join(_pipeline_dir, "validate_results.py"),
                results.output_dir]
        logging.debug("\targs: %s" % (' '.join(map(str, args))))
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode()) 
    #
    # cleanup intermediate files
    #
    msg = "Cleaning up tmp files"
    if keep_tmp:
        skipmsg = "%s because '--keep-tmp' flag set" % msg
        logging.debug("[SKIPPED] %s" % skipmsg)
        shell_commands.append(bash_log(skipmsg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        command = "rm -rf %s" % (results.tmp_dir)
        shell_commands.append(command)    
        shell_commands.append(bash_check_retcode(msg))   
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
            job_id = submit_pbs(shell_commands)
            print job_id
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
    parser.add_argument("--keep-tmp", dest="keep_tmp", action="store_true", 
                        default=False)
    parser.add_argument("library_xml_file")
    parser.add_argument("config_xml_file")
    parser.add_argument("server_name")
    args = parser.parse_args()
    return run(args.library_xml_file, args.config_xml_file, 
               args.server_name, args.num_processors,
               args.stdout_file, args.stderr_file,
               args.dryrun, args.keep_tmp)

if __name__ == '__main__': 
    sys.exit(main())
