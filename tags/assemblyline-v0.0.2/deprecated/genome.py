'''
Created on Mar 11, 2011

@author: mkiyer
'''

def get_genome(genome_version):
    if genome_version == 'hg18':
        return HG18()

class Genome(object):
    pass

class HG18(Genome):
    id = 'hg18'
    chrom_sizes = {'chr1': 247249719,
                   'chr10': 135374737,
                   'chr11': 134452384,
                   'chr12': 132349534,
                   'chr13': 114142980,
                   'chr14': 106368585,
                   'chr15': 100338915,
                   'chr16': 88827254,
                   'chr17': 78774742,
                   'chr18': 76117153,
                   'chr19': 63811651,
                   'chr2': 242951149,
                   'chr20': 62435964,
                   'chr21': 46944323,
                   'chr22': 49691432,
                   'chr3': 199501827,
                   'chr4': 191273063,
                   'chr5': 180857866,
                   'chr6': 170899992,
                   'chr7': 158821424,
                   'chr8': 146274826,
                   'chr9': 140273252,
                   'chrX': 154913754,
                   'chrY': 57772954,
                   'chrM': 16571
                   }        

    def get_chrom_names(self):
        return self.chrom_sizes.keys()

    def get_chrom_sizes(self):
        return self.chrom_sizes

    def get_chrom_size(self, chrom):
        return self.chrom_sizes[chrom]