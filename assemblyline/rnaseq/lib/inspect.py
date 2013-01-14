'''
Created on Aug 10, 2011

@author: mkiyer
'''
import array
import logging
import random

# local imports
import matplotlib
matplotlib.use('Agg')

def inspect_sr_sam(samfh, max_frags=None):
    # process input sam stream
    frags = 0
    unmapped = 0
    useable = 0
    for r in samfh:
        frags += 1
        if r.is_unmapped:
            unmapped += 1
            continue
        useable += 1
        yield -1, r.is_reverse
        if (max_frags is not None) and (useable >= max_frags):
            break
    logging.debug("Processed fragments: %d" % (frags))
    logging.debug("Unmapped fragments: %d" % (unmapped))
    logging.debug("Useable fragments: %d" % (useable))

def inspect_pe_sam(samfh, max_frags=None):
    # function to parse concordant read pairs    
    def parse_pe(infh):
        try:
            while True:
                a = infh.next()
                b = infh.next()
                if a.is_read1:
                    assert b.is_read2
                    r1 = a; r2 = b
                elif a.is_read2:
                    assert b.is_read1
                    r1 = b; r2 = a
                assert r1.qname == r2.qname
                assert r1.is_read1
                assert r2.is_read2
                yield r1,r2
        except StopIteration:
            pass
    # process input sam stream
    frags = 0
    unmapped = 0
    useable = 0
    for r1,r2 in parse_pe(samfh):
        frags += 1
        if r1.is_unmapped or r2.is_unmapped:
            unmapped += 1
            continue
        useable += 1
        # get insert size
        tlen = r1.tlen
        if tlen < 0: tlen = -tlen        
        # yield read info
        yield tlen, r1.is_reverse
        if (max_frags is not None) and (useable >= max_frags):
            break
    logging.debug("Processed fragments: %d" % (frags))
    logging.debug("Unmapped fragments: %d" % (unmapped))
    logging.debug("Useable fragments: %d" % (useable))

class RnaseqLibraryMetrics(object):
    def __init__(self, min_frag_size=0, max_frag_size=0, cutoff_frac=0.90):
        self.min_frag_size = max(0, min_frag_size)
        self.max_frag_size = max(0, max_frag_size)
        self.arr = array.array('L', (0 for x in xrange(self.min_frag_size, self.max_frag_size+1)))
        self.read1_rev_count = 0
        self.read1_count = 0
        self.cutoff_frac = min(1.0, max(0.0, cutoff_frac))

    def clear_frag_size(self):
        self.arr = array.array('L', (0 for x in xrange(self.min_frag_size, self.max_frag_size+1)))
    def clear_strandedness(self):
        self.read1_count = 0
        self.read1_rev_count = 0

    def from_stream(self, myiter):
        total = 0
        inrange = 0
        for tlen, is_reverse in myiter:
            if tlen >= 0:
                if (self.min_frag_size <= tlen <= self.max_frag_size):
                    self.arr[tlen - self.min_frag_size] += 1
                    inrange += 1
            self.read1_rev_count += int(is_reverse)
            self.read1_count += 1
            total += 1
        logging.debug("Stream had %d total samples with %d within the fragment size range" % (total, inrange))

    def frag_size_from_random(self, mean, stdev, max_frags=100000):
        total = 0
        inrange = 0
        while total < max_frags:
            tlen = int(round(random.normalvariate(mean, stdev),0))
            if (self.min_frag_size <= tlen <= self.max_frag_size):
                self.arr[tlen - self.min_frag_size] += 1
                inrange += 1
            total += 1
        logging.debug("Randomly sampled %d insert sizes with %d within the fragment size range" % (total, inrange))

    def strandedness_from_random(self, max_frags=100000):
        n = 0
        nrev = 0
        while n < max_frags:
            rev = random.randint(0,1)
            nrev += rev
            n += 1
            self.read1_rev_count += rev
            self.read1_count += 1
        logging.debug("Randomly sampled %d orientations with %d reverse" % (n,nrev))

    def read1_strand_fraction(self):
        return self.read1_rev_count / float(self.read1_count)

    def predict_library_type(self):
        frac = self.read1_rev_count / float(self.read1_count)
        if frac >= self.cutoff_frac:
            return 'fr-firststrand'
        elif frac <= (1.0 - self.cutoff_frac):
            return 'fr-secondstrand'
        else:
            return 'fr-unstranded'

    def tlen_at_percentile(self, per):
        n = sum(self.arr)
        per_n = n * per / 100.0
        count = 0
        for tlen,x in enumerate(self.arr): 
            count += x
            if (count >= per_n):
                break
        return tlen + self.min_frag_size
    
    def percentile_at_tlen(self, tlen):
        if tlen < self.min_frag_size:
            return 0.0
        elif tlen > self.max_frag_size:
            return 100.0
        ind = tlen - self.min_frag_size
        count_le = sum(self.arr[:ind+1])        
        per = 100.0 * count_le / float(sum(self.arr))
        return per

    def num_frag_size_samples(self):
        return sum(self.arr)
    
    def mode(self):
        return self.arr.index(max(self.arr)) + self.min_frag_size

    def mean(self):
        count = 0
        n = 0        
        for i,x in enumerate(self.arr): 
            count += i*x
            n += x
        if n == 0:
            return None            
        return self.min_frag_size + (count / float(n))
   
    def std(self):
        # first get mean (not shifted)
        mean = self.mean()
        if mean is None:
            return None
        n = 0
        std = 0
        for i,x in enumerate(self.arr):
            tlen = i + self.min_frag_size
            std = std + x*((tlen - mean)**2)
            n += x
        std = (std / float(n-1))**0.5
        return std

    def plot_frag_size_dist(self, filename, max_percentile=100.0):
        import matplotlib.pyplot as plt
        plot_max_tlen = self.tlen_at_percentile(max_percentile)
        xdata = range(self.min_frag_size, plot_max_tlen+1)
        num_samples = self.num_frag_size_samples()
        ydata = [float(self.arr[x-self.min_frag_size])/num_samples for x in xdata] 
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata,ydata)
        ax.minorticks_on()
        plt.title("Fragment size distribution")
        plt.xlabel("Fragment size (bp)")
        plt.ylabel("Fraction of reads")
        plt.savefig(filename)
        plt.close()

    def to_file(self, fileh):
        print >>fileh, '\t'.join(["#read1_rev_count", "read1_count"])
        print >>fileh, '\t'.join(map(str, [self.read1_rev_count, self.read1_count]))
        print >>fileh, '\t'.join(["#frac", "cutoff_frac", "library_type"])
        print >>fileh, '\t'.join(map(str, [self.read1_strand_fraction(),
                                           self.cutoff_frac,
                                           self.predict_library_type()]))
        print >>fileh, '\t'.join(["#insert_size", "num_samples"])
        for i,x in enumerate(self.arr):
            print >>fileh, '\t'.join([str(i + self.min_frag_size), str(x)])        

    @staticmethod
    def from_file(filename):
        res = RnaseqLibraryMetrics()
        # skip first comment
        fileh = open(filename)
        fileh.next()
        # read strandedness
        fields = map(int, fileh.next().strip().split('\t'))
        res.read1_rev_count = fields[0]
        res.read1_count = fields[1]
        # skip comment
        fileh.next()
        # read cutoff frac
        fields = fileh.next().strip().split('\t')
        res.cutoff_frac = float(fields[1])
        # read frag size distribution
        tlens = []
        counts = []
        for line in fileh:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            i,x = map(int, fields[0:2])
            tlens.append(i)
            counts.append(x)
        res.min_frag_size = tlens[0]
        res.max_frag_size = tlens[-1]
        res.arr = array.array('L', counts)
        fileh.close()
        return res
