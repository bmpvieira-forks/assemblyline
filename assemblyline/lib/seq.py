'''
Created on Nov 1, 2011

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
import string

#Translation table for reverse Complement, with ambiguity codes
DNA_COMPLEMENT = string.maketrans( "ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb" )
RNA_COMPLEMENT = string.maketrans( "ACGURYKMBDHVacgurykmbdhv", "UGCAYRMKVHDBugcayrmkvhdb" )
#Translation table for DNA <--> RNA
DNA_TO_RNA = string.maketrans( "Tt", "Uu" )
RNA_TO_DNA = string.maketrans( "Uu", "Tt" )

def DNA_complement( sequence ):
    '''complement DNA sequence string'''
    return sequence.translate( DNA_COMPLEMENT )
def DNA_reverse_complement( sequence ):
    '''returns the reverse complement of the sequence'''
    return DNA_complement(sequence[::-1])
def to_DNA( sequence ):
    return sequence.translate( DNA_TO_RNA )
#complement RNA sequence string
def RNA_complement( sequence ):
    return sequence.translate( RNA_COMPLEMENT )
def RNA_reverse_complement( self, sequence ):
    return RNA_complement( sequence[::-1] )
def to_RNA( sequence ):
    return sequence.translate( RNA_TO_DNA )