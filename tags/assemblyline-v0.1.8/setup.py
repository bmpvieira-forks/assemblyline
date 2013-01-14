'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012,2013 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# ---- Extension Modules ----------------------------------------------------
def get_extension_modules():
    extensions = []
    # Interval clustering 
    extensions.append(Extension("assemblyline.lib.bx.cluster", 
                                ["assemblyline/lib/bx/cluster.pyx", 
                                 "assemblyline/lib/bx/intervalcluster.c"], 
                                include_dirs=["assemblyline/lib/bx"]))
    # Interval intersection
    extensions.append(Extension("assemblyline.lib.bx.intersection", 
                                ["assemblyline/lib/bx/intersection.pyx"]))
    return extensions

def main():
    setup(name="assemblyline",
          version='0.1.6',
          description="transcriptome meta-assembly from RNA-Seq",
          long_description=__doc__,
          author = "Matthew Iyer",
          author_email = "mkiyer@umich.edu",
          license="GPL3",
          platforms="Linux",
          url="http://assemblyline.googlecode.com",
          packages=['assemblyline',
                    'assemblyline.lib',
                    'assemblyline.lib.bx',
                    'assemblyline.lib.assemble',
                    'assemblyline.pipeline',
                    'assemblyline.test',
                    'assemblyline.utils',
                    'assemblyline.rnaseq',
                    'assemblyline.rnaseq.lib',
                    'assemblyline.rnaseq.pipeline',
                    'assemblyline.rnaseq.utils'],
          package_data={'assemblyline.lib': ['classify_transcripts.R']},                             
          ext_modules=get_extension_modules(),
          cmdclass= {'build_ext': build_ext})

if __name__ == '__main__': 
    main()