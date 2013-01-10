'''
Created on Dec 17, 2011

@author: mkiyer

AssemblyLine: transcriptome meta-assembly from RNA-Seq

Copyright (C) 2012 Matthew Iyer

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
from assemblyline.lib.transcript import Exon

# transcript graph attributes
TRANSCRIPT_IDS = 'tids'
STRAND_SCORE = 'strand_score'

# assembler graph attributes
NODE_SCORE = 'score'
NODE_LENGTH = 'length'
SMOOTH_FWD = 'smfwd'
SMOOTH_REV = 'smrev'
SMOOTH_TMP = 'smtmp'
CHAIN_NODES = 'chain'
CHAIN_DATA = 'chain_data'
CHAIN_EDGES = 'chain_edges'

# constant minimum path score
MIN_SCORE = 1.0e-8


