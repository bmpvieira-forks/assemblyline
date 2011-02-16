'''
Created on Nov 29, 2010

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# ---- Extension Modules ----------------------------------------------------
def get_extension_modules():
    extensions = []
    # Graph node
    #extensions.append( Extension( "lib.cnode", ["lib/cnode.pyx"] ) )
    # Interval clustering                
    extensions.append( Extension( "lib.bx.cluster", [ "lib/bx/cluster.pyx", "lib/bx/intervalcluster.c"], 
                                  include_dirs=["lib/bx"]) )
    # Interval intersection
    extensions.append( Extension( "lib.bx.intersection", [ "lib/bx/intersection.pyx" ] ) )
    return extensions

def main():
    setup(name = "assemblyline",
          ext_modules = get_extension_modules(),
          author = "Matthew Iyer",
          author_email = "mkiyer@umich.edu",
          description = "Transcriptome assembly clustering tool",
          url = "http://assemblyline.googlecode.com",
          cmdclass= {'build_ext': build_ext})

if __name__ == '__main__': main()