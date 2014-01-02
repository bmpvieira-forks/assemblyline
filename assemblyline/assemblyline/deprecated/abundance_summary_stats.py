'''
Created on Sep 18, 2013

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import numpy as np

class ExpressionMatrix(object):
    def __init__(self):
        self.pheno_header_fields = None
        self.pheno_lines = []
        self.meta_header_fields = None
        self.meta_lines = []
        self.mat = None

    @staticmethod
    def open(matrix_dir):
        m = ExpressionMatrix()
        pheno_file = os.path.join(matrix_dir, 'phenos.txt')
        metadata_file = os.path.join(matrix_dir, 'metadata.txt')
        matrix_file = os.path.join(matrix_dir, 'isoform_fpkm.mmap')
        # read pheno file
        fileh = open(pheno_file)
        m.pheno_header_fields = fileh.next().strip().split('\t')
        m.pheno_lines = []
        for line in fileh:
            fields = line.strip().split('\t')
            m.pheno_lines.append(fields)
        fileh.close()
        # read metadata file
        fileh = open(metadata_file)
        m.meta_header_fields = fileh.next().strip().split('\t')
        m.meta_lines = []
        for line in fileh:
            fields = line.strip().split('\t')
            m.meta_lines.append(fields)
        fileh.close()
        # read matrix data
        m.mat = np.memmap(matrix_file, dtype='float32', mode='r', 
                          shape=(len(m.meta_lines),
                                 len(m.pheno_lines)))
        return m

    def close(self):
        pass

    def get_transcript_expression(self):
        pass


def main():
    # Command line parsing
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('matrix_dir')
    args = parser.parse_args()
    if not os.path.exists(args.matrix_dir):
        parser.error("Expression matrix directory '%s' not found" % (args.matrix_dir))
    # read matrix
    m = ExpressionMatrix.open(args.matrix_dir)
    #category_ind = m.meta_header_fields.index['category']
    header_fields = list(m.meta_header_fields)
    header_fields.extend(['avg', 'q0', 'q1', 'q5', 'q10', 'q25', 'q50', 
                          'q75', 'q90', 'q95', 'q99', 'q100'])
    print '\t'.join(header_fields)
    q = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]
    for i in xrange(m.mat.shape[0]):
        a = m.mat[i,:]
        fields = list(m.meta_lines[i])
        fields.append(np.mean(a))
        fields.extend(np.percentile(a, q))
        print '\t'.join(map(str, fields))
        
    return 0


if __name__ == '__main__':
    sys.exit(main())
