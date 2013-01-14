'''
Created on Oct 11, 2012

@author: mkiyer
'''
import logging
import argparse
import sys

import assemblyline.rnaseq.lib.config as config

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("config_xml_file")
    parser.add_argument("server_name")
    args = parser.parse_args()
    # read and validate configuration file
    logging.info("Reading pipeline configuration file")
    pipeline = config.PipelineConfig.from_xml(args.config_xml_file)
    if not pipeline.is_valid(args.server_name):
        logging.error("Pipeline config not valid")
        return config.JOB_ERROR
    if not pipeline.is_software_valid():
        logging.error("Pipeline software not valid")
        return config.JOB_ERROR
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())