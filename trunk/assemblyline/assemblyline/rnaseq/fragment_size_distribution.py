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

# project imports
from sam import parse_pe_reads

class FragmentSizeDistribution(object):
    def __init__(self):
        self.min_isize = 0
        self.max_isize = 0
        self.arr = array.array('L', (0 for x in xrange(self.min_isize, self.max_isize+1)))

    def isize_at_percentile(self, per):
        n = sum(self.arr)
        per_n = n * per / 100.0
        count = 0
        for isize,x in enumerate(self.arr): 
            count += x
            if (count >= per_n):
                break
        return isize + self.min_isize
    
    def percentile_at_isize(self, isize):
        if isize < self.min_isize:
            return 0.0
        elif isize > self.max_isize:
            return 100.0
        ind = isize - self.min_isize
        count_le = sum(self.arr[:ind+1])        
        per = 100.0 * count_le / float(sum(self.arr))
        return per

    @property
    def n(self):
        if self.arr is None: return 0
        return sum(self.arr)
    
    def mode(self):
        return self.arr.index(max(self.arr)) + self.min_isize

    def mean(self):
        count = 0
        n = 0        
        for i,x in enumerate(self.arr): 
            count += i*x
            n += x
        if n == 0:
            return None            
        return self.min_isize + (count / float(n))
   
    def std(self):
        # first get mean (not shifted)
        mean = self.mean()
        if mean is None:
            return None
        n = 0
        std = 0
        for i,x in enumerate(self.arr):
            isize = i + self.min_isize
            std = std + x*((isize - mean)**2)
            n += x
        std = (std / float(n-1))**0.5
        return std
    
    @staticmethod
    def merge(a, b):
        """
        merge two fragment size distributions
        """
        min_isize = min(a.min_isize, b.min_isize)
        max_isize = max(a.max_isize, b.max_isize)
        arr = array.array('L', (0 for x in xrange(min_isize, max_isize+1)))
        # add counts from both distributions
        for isize in xrange(min_isize, max_isize):
            count = 0
            # add from this array
            if a.min_isize <= isize <= a.max_isize:
                count += a.arr[isize - a.min_isize]
            # add from b array
            if b.min_isize <= isize <= b.max_isize:
                count += b.arr[isize - b.min_isize]
            # store count
            arr[isize - min_isize] = count
        d = FragmentSizeDistribution()
        d.min_isize = min_isize
        d.max_isize = max_isize
        d.arr = arr
        return d

    def plot(self, filename, max_percentile=100.0):
        import matplotlib.pyplot as plt
        plot_max_isize = self.isize_at_percentile(max_percentile)
        xdata = range(self.min_isize, plot_max_isize+1)
        num_samples = self.n
        ydata = [float(self.arr[x-self.min_isize])/num_samples for x in xdata] 
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
        print >>fileh, '\t'.join(["#insert_size", "num_samples"])
        for i,x in enumerate(self.arr):
            print >>fileh, '\t'.join([str(i + self.min_isize), str(x)])        

    @staticmethod
    def from_file(fileh):
        isizes = []
        counts = []
        for line in fileh:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            i,x = map(int, fields[0:2])
            isizes.append(i)
            counts.append(x)
        res = FragmentSizeDistribution()
        res.min_isize = isizes[0]
        res.max_isize = isizes[-1]
        res.arr = array.array('L', counts)
        return res 

    @staticmethod
    def from_random(mean, stdev, min_isize, max_isize, samples=100000):
        """
        initialize from a random sample using normal distribution with 
        mean 'mean' and stdev 'stdev'
        """
        d = FragmentSizeDistribution()
        # implement simple checks
        assert min_isize < mean < max_isize
        assert stdev < (max_isize - min_isize)
        # initialize
        d.min_isize = min_isize
        d.max_isize = max_isize
        d.arr = array.array('L', (0 for x in xrange(min_isize, max_isize+1)))
        count = 0
        outside_range = 0
        while True:
            if count > samples:
                break
            isize = int(round(random.normalvariate(mean, stdev),0))
            if (min_isize <= isize <= max_isize):
                # store in array
                d.arr[isize - min_isize] += 1
                count += 1
            else:
                outside_range += 1
        return d

    @staticmethod
    def from_sam(samfh, min_isize, max_isize, max_samples=None):
        """
        iterates through a SAM/BAM object looking for uniquely mapping 
        concordant reads.  keeps a histogram of all observed insert 
        sizes in the reads.  stops once 'max_samples' valid reads are 
        encountered, or the end of the file is reached
        """
        d = FragmentSizeDistribution()
        d.min_isize = min_isize
        d.max_isize = max_isize
        d.arr = array.array('L', (0 for x in xrange(min_isize, max_isize+1)))        
        frags = 0
        good_frags = 0
        unmapped = 0
        ambiguous = 0
        outside_range = 0
        for pe_reads in parse_pe_reads(samfh):
            frags += 1
            # only allow mappings where there is a single
            # insert size (multiple isoforms are ambiguous)
            isizes = set()        
            for r in pe_reads[0]:
                if r.is_unmapped:
                    continue
                # get insert size
                isize = r.isize
                if isize < 0: isize = -isize
                isizes.add(isize)
            # insert size must be within range
            if len(isizes) == 0:
                unmapped += 1
            elif len(isizes) > 1:
                ambiguous += 1
            else:
                isize = isizes.pop()
                if (min_isize <= isize <= max_isize):
                    # store in array
                    d.arr[isize - min_isize] += 1
                    good_frags += 1
                else:
                    outside_range += 1
            if (max_samples is not None) and (good_frags >= max_samples):
                break
        logging.debug("Processed fragments: %d" % (frags))
        logging.debug("Unambiguous paired fragments: %d" % (good_frags))
        logging.debug("Unmapped fragments: %d" % (unmapped))
        logging.debug("Ambiguous fragments: %d" % (ambiguous))
        logging.debug("Fragments outside range: %d" % (outside_range))
        return d

