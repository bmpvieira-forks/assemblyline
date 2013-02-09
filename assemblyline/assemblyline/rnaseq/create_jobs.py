'''
Created on Aug 7, 2011

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import subprocess
import shutil

import assemblyline.rnaseq.lib.config as config
from assemblyline.rnaseq.lib.base import detect_format
from assemblyline.rnaseq.lib.libtable import Library, read_library_table_xls, FRAGMENT_LAYOUT_PAIRED
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
    commands.append("module purge")
    for name in pipeline.modules:
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

def create_job(library, pipeline, server, config_xml_file, 
               num_processors, keep_tmp, overwrite):
    # search for library sequence files
    success = config.resolve_library_sequence_files(server, library)
    if not success:
        logging.error("Library %s sequence files not found" % (library.library_id))
        return config.JOB_ERROR
    # check for valid data structure
    if not library.is_valid():
        logging.error("Library %s not valid" % (library.library_id))
        return config.JOB_ERROR
    # check genome
    if library.species not in pipeline.genomes:
        logging.error("Library %s species %s not supported" % 
                      (library.library_id, library.species))
        return config.JOB_ERROR
    # resolve genome paths now that server is known
    genome_local = pipeline.genomes[library.species]
    genome_static = genome_local.resolve_paths(server.references_dir)    
    # setup results
    output_dir = os.path.join(server.output_dir, library.library_id)
    results = config.RnaseqResults(library, pipeline, output_dir)
    # create output directories
    if os.path.exists(results.output_dir):
        if not overwrite:
            logging.error("Library %s output dir %s exists and --overwrite "
                          "option not set" % 
                          (library.library_id, results.output_dir))
            return config.JOB_ERROR
    else:
        logging.debug("Creating directory: %s" % (results.output_dir))
        os.makedirs(results.output_dir)
    if not os.path.exists(results.tmp_dir):
        os.makedirs(results.tmp_dir)
    if not os.path.exists(results.log_dir):
        os.makedirs(results.log_dir)
    #
    # format library information as xml
    #
    logging.debug("Writing library information to XML")
    library.to_xml_file(results.library_xml_file)
    #
    # copy configuration xml file
    #
    logging.debug("Copying configuration file")
    shutil.copyfile(config_xml_file, results.config_xml_file)
    #
    # build up sequence of commands
    #
    shell_commands = []
    #
    # tmp directory may have been removed at end of a previous run
    #
    if not os.path.exists(results.tmp_dir):
        logging.debug("Creating directory: %s" % (results.tmp_dir))
        shell_commands.append("mkdir -p %s" % (results.tmp_dir))
        shell_commands.append(bash_check_retcode())   
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
            results.config_xml_file, server.name]
    command = ' '.join(map(str, args))
    logging.debug("\tcommand: %s" % (command))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))   
    #
    # extract fastq from bam input files
    #
    msg = "Converting BAM input files to FASTQ"
    if len(results.bam_fastq_prefixes) > 0:    
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'extract_fastq_from_bam.log')
        # convert bam to fastq
        for i,prefix in enumerate(results.bam_fastq_prefixes):
            args = ["python",
                    os.path.join(_pipeline_dir, "extract_fastq_from_bam.py"),
                    "--no-unpaired",
                    "--readnum-in-qname",
                    "--tmp-dir", results.tmp_dir,
                    library.bam_files[i], prefix,
                    '> %s 2>&1' % (log_file)] 
            command = ' '.join(map(str, args))
            shell_commands.append(command)
            shell_commands.append(bash_check_retcode())
    #
    # copy/concatenate sequences read1
    #
    msg = "Concatenating/copying read1 sequence files"
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    commands = concatenate_sequences(library.read1_files + results.bam_read1_files, 
                                     results.copied_fastq_files[0])
    shell_commands.extend(commands)
    #
    # copy/concatenate sequences read2
    #
    if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:    
        msg = "Concatenating/copying read2 sequence files"
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        commands = concatenate_sequences(library.read2_files + results.bam_read2_files, 
                                         results.copied_fastq_files[1])
        shell_commands.extend(commands)
    #
    # Run FASTQC quality assessment
    #
    msg = "Running FASTQC quality assessment"
    logging.debug("%s" % msg)    
    shell_commands.append(bash_log(msg, "INFO"))
    num_threads = min(num_processors, len(results.copied_fastq_files))
    log_file = os.path.join(results.log_dir, 'fastqc.log')
    args = ['fastqc', "--threads", num_threads, "-o", results.output_dir]
    args.extend(results.copied_fastq_files)
    args.append('> %s 2>&1' % (log_file))
    logging.debug("\targs: %s" % (' '.join(map(str, args))))
    command = ' '.join(map(str, args))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode())   
    #
    # Filter reads that map to abundant sequences
    #
    msg = "Filtering reads that map to abundant sequences"
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    args = ["python", 
            os.path.join(_pipeline_dir, "filter_reads.py"),
            "--num-processors", num_processors,
            "--tmp-dir", results.tmp_dir,
            genome_static.abundant_bowtie2_index,
            ','.join(results.fastqc_data_files),
            ','.join(results.copied_fastq_files),
            ','.join(results.filtered_fastq_files),
            results.sorted_abundant_bam_file,
            results.abundant_counts_file]
    logging.debug("\targs: %s" % (' '.join(map(str, args))))
    command = ' '.join(map(str, args))
    log_file = os.path.join(results.log_dir, 'filter_abundant_sequences.log')
    command += ' > %s 2>&1' % (log_file)
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode())
    #
    # inspect the library and determine characteristics
    #
    msg = "Inspecting library"
    logging.debug(msg) 
    shell_commands.append(bash_log(msg, "INFO"))
    args = ["python",
            os.path.join(_pipeline_dir, "inspect_library.py"),
            "-p", num_processors]
    for arg in pipeline.inspect_args:
        args.extend(arg.split())
    args.extend([genome_static.fragment_size_bowtie1_index,
                 results.library_metrics_file,
                 results.frag_size_dist_plot_file])
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
    if (not pipeline.tophat_fusion_run):
        logging.debug("[SKIPPED] %s" % (msg))
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        msg = "Aligning reads with Tophat-Fusion"
        logging.debug("%s" % (msg))
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_tophat.py"),
                "-p", num_processors,
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
        # index tophat fusion bam file
        msg = "Indexing Tophat-Fusion BAM file"
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
    msg = "Post-processing Tophat-Fusion results"
    skip = ((not pipeline.tophat_fusion_run) and
            (not pipeline.tophat_fusion_post_run)) 
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
        shell_commands.extend(bash_remove_files(results.tophat_fusion_post_tmp_files))
    #
    # align reads with tophat
    #
    msg = "Aligning reads with Tophat"
    logging.debug("%s" % (msg))
    shell_commands.append(bash_log(msg, "INFO"))
    args = ["python", os.path.join(_pipeline_dir, "run_tophat.py"),
            "-p", num_processors,
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
    msg = "Indexing Tophat BAM file"
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
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    args = ["java", "-Xmx4g", "-jar", 
            "$PICARDPATH/CollectMultipleMetrics.jar",
            "INPUT=%s" % (results.tophat_bam_file),
            "REFERENCE_SEQUENCE=%s" % genome_static.genome_fasta_file,
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
            genome_static.genome_fasta_file,   
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
    msg = "Converting unmapped BAM to FASTQ for realignment"
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    log_file = os.path.join(results.log_dir, 'unmapped_bam_to_fastq.log')
    # convert bam to paired fastq
    args = ["python",
            os.path.join(_pipeline_dir, "extract_fastq_from_bam.py"),
            "--tmp-dir", results.tmp_dir,
            results.tophat_unmapped_bam_file,
            results.unmapped_fastq_prefix,
            '> %s 2>&1' % (log_file)]
    command = ' '.join(map(str, args))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode(msg))
    #
    # align to pathogen sequences
    #
    msg = "Aligning to pathogen sequences"
    logging.debug(msg) 
    shell_commands.append(bash_log(msg, "INFO"))
    log_file = os.path.join(results.log_dir, 'bowtie2_align_pathogens.log')        
    args = ["python",
            os.path.join(_pipeline_dir, "bowtie2_align.py"),
            "-p", num_processors,
            "--tmp-dir", results.tmp_dir,
            "-1", results.unmapped_paired_fastq_files[0],
            "-2", results.unmapped_paired_fastq_files[1],
            "-U", ','.join(results.unmapped_unpaired_fastq_files)]
    # extra args
    for arg in pipeline.pathogen_screen_bt2_args:
        args.append('--extra-arg="%s"' % arg)
    args.extend([genome_static.pathogen_bowtie2_index,                     
                 results.pathogen_bam_file])
    args.append('> %s 2>&1' % (log_file))
    command = ' '.join(map(str, args))
    logging.debug("\tcommand: %s" % (command))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode())
    #
    # count pathogen sequences
    #
    msg = "Counting pathogen alignments"
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    log_file = os.path.join(results.log_dir, 'pathogen_counts.log')
    args = ["python",
            os.path.join(_pipeline_dir, "bowtie2_quantify.py"),
            results.pathogen_bam_file,
            results.pathogen_counts_file,
            '> %s 2>&1' % (log_file)]
    command = ' '.join(map(str, args))
    logging.debug("\tcommand: %s" % (command))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode())            
    #
    # align to repeat elements
    #
    msg = "Aligning to repeat element sequences"
    logging.debug(msg) 
    shell_commands.append(bash_log(msg, "INFO"))
    log_file = os.path.join(results.log_dir, 'bowtie2_align_repeats.log')        
    args = ["python",
            os.path.join(_pipeline_dir, "bowtie2_align.py"),
            "-p", num_processors,
            "--tmp-dir", results.tmp_dir,
            "-1", results.unmapped_paired_fastq_files[0],
            "-2", results.unmapped_paired_fastq_files[1],
            "-U", ','.join(results.unmapped_unpaired_fastq_files)]
    # extra args
    for arg in pipeline.repeat_elements_bt2_args:
        args.append('--extra-arg="%s"' % arg)
    args.extend([genome_static.repbase_bowtie2_index,                     
                 results.repeat_element_bam_file])
    args.append('> %s 2>&1' % (log_file))
    command = ' '.join(map(str, args))
    logging.debug("\tcommand: %s" % (command))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode())
    #
    # count repeat elements
    #
    msg = "Counting repeat element alignments"
    logging.debug(msg)
    shell_commands.append(bash_log(msg, "INFO"))
    log_file = os.path.join(results.log_dir, 'repeat_element_counts.log')
    args = ["python",
            os.path.join(_pipeline_dir, "bowtie2_quantify.py"),
            results.repeat_element_bam_file,
            results.repeat_element_counts_file,
            '> %s 2>&1' % (log_file)]
    command = ' '.join(map(str, args))
    logging.debug("\tcommand: %s" % (command))
    shell_commands.append(command)
    shell_commands.append(bash_check_retcode()) 
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
    if (not pipeline.cufflinks_ab_initio_run):
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_cufflinks.py"),
                "-p", num_processors,
                "-L", library.library_id]
        if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
            args.append("--learn-frag-size")
        # resolve genome-specific args
        for arg in pipeline.cufflinks_ab_initio_args:
            args.append('--cufflinks-arg="%s"' % genome_static.resolve_arg(arg))
        args.extend([results.tophat_bam_file,
                     results.library_metrics_file,
                     results.cufflinks_ab_initio_dir])
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
    if (not pipeline.cufflinks_known_run):
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_cufflinks.py"),
                "-p", num_processors,
                "-L", library.library_id]
        if library.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
            args.append("--learn-frag-size")
        # resolve genome-specific args
        for arg in pipeline.cufflinks_known_args:
            args.append('--cufflinks-arg="%s"' % genome_static.resolve_arg(arg))
        args.extend([results.tophat_bam_file,
                     results.library_metrics_file,                     
                     results.cufflinks_known_dir])
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
    if (not pipeline.htseq_count_run):
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
    if (not pipeline.varscan_run):
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
    if (not pipeline.varscan_run):
        logging.debug("[SKIPPED] %s" % msg)
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        logging.debug(msg)
        shell_commands.append(bash_log(msg, "INFO"))
        args = ["python", os.path.join(_pipeline_dir, "run_varscan.py")]
        for arg in pipeline.varscan_args:
            args.append('--varscan-arg="%s"' % (arg))
        args.extend(["$VARSCANPATH/VarScan.jar", "mpileup2snp",
                     genome_static.genome_fasta_file,
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
    if (not pipeline.varscan_run):
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
    # convert VCF to annovar input format
    #
    msg = "Converting VCF to Annovar format"
    skip = ((not pipeline.varscan_run) and 
            (not pipeline.annovar_var))
    if skip:
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'convert2annovar.log')
        args = ["perl", "$ANNOVARPATH/convert2annovar.pl",
                "-format", "vcf4", "-filter", "pass", 
                "--includeinfo", "--withzyg", 
                results.varscan_snv_file,
                '> %s 2> %s' % (results.annovar_input_file, log_file)]
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())    
    #
    # run annovar to filter variant calls
    #
    msg = "Annotating variants with Annovar"
    skip = ((not pipeline.varscan_run) and 
            (not pipeline.annovar_var))
    if skip:
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'annovar_summary.log')
        args = ["perl $ANNOVARPATH/summarize_annovar.pl",
                "--buildver", genome_static.annovar_buildver,
                genome_static.annovar_summary_args,
                "--outfile", results.annovar_output_prefix,
                results.annovar_input_file,
                genome_static.annovar_db,
                '> %s 2>&1' % (log_file)]
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # annotate cosmic positions
    #
    msg = "Annotating COSMIC positions with Annovar"
    skip = ((not pipeline.varscan_run) and 
            (not pipeline.annovar_var))
    if skip:
        shell_commands.append(bash_log(msg, "SKIPPED"))
    else:
        shell_commands.append(bash_log(msg, "INFO"))
        log_file = os.path.join(results.log_dir, 'annovar_cosmic.log')
        args = ["perl $ANNOVARPATH/annotate_variation.pl", 
                "--buildver", genome_static.annovar_buildver,
                "--filter",
                "--dbtype", genome_static.annovar_cosmicver,
                results.annovar_input_file,
                genome_static.annovar_db,
                '> %s 2>&1' % (log_file)]
        command = ' '.join(map(str, args))
        shell_commands.append(command)
        shell_commands.append(bash_check_retcode())
    #
    # write job finished file
    #
    msg = "Validating results"
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
    # write shell/pbs scripts
    #
    if server.pbs:
        pbsfh = open(results.pbs_script_file, "w")
        print >>pbsfh, "#!/bin/sh"
        pbs_commands = get_pbs_header(job_name=library.library_id,
                                      num_processors=num_processors,
                                      node_processors=server.node_processors,
                                      node_memory=server.node_mem,
                                      pbs_script_lines=server.pbs_script_lines, 
                                      walltime=config.PBS_JOB_WALLTIME, 
                                      mem=config.PBS_JOB_MEM,
                                      stdout_filename=results.pbs_stdout_file,
                                      stderr_filename=results.pbs_stderr_file)    
        for command in pbs_commands:
            print >>pbsfh, command
        for command in shell_commands:
            print >>pbsfh, command
        pbsfh.close()
    shellfh = open(results.shell_script_file, "w")
    print >>shellfh, "#!/bin/sh"
    for command in shell_commands:
        print >>shellfh, command
    shellfh.close()
    return config.JOB_SUCCESS

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", 
                        dest="verbose", default=False,
                        help="show debugging logging messages "
                        "[default=%(default)s]")
    parser.add_argument("-p", type=int, dest="num_processors", 
                        default=1,
                        help="Number of processors per job")
    parser.add_argument("--keep-tmp", dest="keep_tmp", 
                        action="store_true", default=False,
                        help="Retain intermediate files generated during "
                        "analysis (warning: may consume large amounts of "
                        "disk space, only recommended for debugging "
                        "purposes")
    parser.add_argument("-f", "--overwrite", action="store_true",
                        default=False,
                        help="overwrite existing job files "
                        "[default=%(default)s]")
    inp_group = parser.add_mutually_exclusive_group(required=True)
    inp_group.add_argument("--xls", dest="library_xls_file", 
                           default=None,
                           help="Excel (.xls/.xlsx) spreadsheet "
                           "containing RNA-Seq library information")
    inp_group.add_argument("--xml", dest="library_xml_file", 
                           default=None,
                           help="XML formatted file containing RNA-Seq "
                           "library information")
    parser.add_argument("--library-id", dest="input_library_ids",
                        action="append", default=None,
                        help="specific library identifier(s) "
                        "to process [entire XML/Excel files processed "
                        "by default]")
    parser.add_argument("config_xml_file")
    parser.add_argument("server_name")
    args = parser.parse_args()
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("Assemblyline version %s" % (assemblyline.__version__))
    logging.info("=============================")
    #
    # read and validate configuration file
    #
    if not os.path.exists(args.config_xml_file):
        parser.error("Configuration XML file %s not found" % 
                     (args.config_xml_file))
    logging.info("Reading pipeline configuration file")
    pipeline = config.PipelineConfig.from_xml(args.config_xml_file)
    if not pipeline.is_valid(args.server_name):
        logging.error("Pipeline config not valid")
        return config.JOB_ERROR
    server = pipeline.servers[args.server_name]
    # sublist of library ids
    if args.input_library_ids is None:
        input_library_ids = None
    else:
        input_library_ids = set(args.input_library_ids)
    #
    # read libraries
    #
    if args.library_xml_file is not None:
        if not os.path.exists(args.library_xml_file):
            parser.error("XML file %s not found" % (args.library_xml_file))
        # read library from XML
        libraries = {}
        for library in Library.from_xml_file(args.library_xml_file):
            if ((input_library_ids is None) or 
                (library.library_id in input_library_ids)):
                libraries[library.library_id] = library
    elif args.library_xls_file is not None:                
        if not os.path.exists(args.library_xls_file):
            parser.error("Excel file %s not found" % (args.library_xls_file))
        # read libraries from XLS/XLSX
        logging.info("Reading library table '%s'" % (args.library_xls_file))
        libraries = read_library_table_xls(args.library_xls_file)
        if input_library_ids is not None:            
            libraries = dict((x,libraries[x]) for x in input_library_ids)
    #
    # process each library
    #
    logging.info("Generating jobs")
    for library in libraries.itervalues():
        logging.info("Library: %s" % (library.library_id)) 
        retcode = create_job(library, pipeline, server, 
                             args.config_xml_file, 
                             args.num_processors, 
                             args.keep_tmp,
                             args.overwrite)
        if retcode != 0:
            logging.error("Library %s error" % (library.library_id))
    logging.info("Done")
    return config.JOB_SUCCESS

if __name__ == '__main__': 
    sys.exit(main())
