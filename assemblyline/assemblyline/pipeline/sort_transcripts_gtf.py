'''
Created on Nov 4, 2010

@author: mkiyer
'''
import sys
from assemblyline.lib.gtf import sort_gtf

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    sort_gtf(input_file, output_file)    
