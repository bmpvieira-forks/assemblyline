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
    for r in samfh:
        frags += 1
        if r.is_unmapped:
            unmapped += 1
            continue
        yield None, r.is_reverse
        if (max_frags is not None) and (frags >= max_frags):
            break
    logging.debug("Processed fragments: %d" % (frags))
    logging.debug("Unmapped fragments: %d" % (unmapped))

def inspect_pe_sam(samfh, max_frags=None):
    # function to parse concordant read pairs    
    def parse_pe(infh):
        try:
            while True:
                r1 = infh.next()
                r2 = infh.next()
                assert r1.qname == r2.qname
                assert r1.is_read1
                assert r2.is_read2
                yield r1,r2
        except StopIteration:
            pass
    # process input sam stream
    frags = 0
    unmapped = 0
    for r1,r2 in parse_pe(samfh):
        frags += 1
        if r1.is_unmapped or r2.is_unmapped:
            unmapped += 1
            continue
        # get insert size
        isize = r1.isize
        if isize < 0: isize = -isize        
        # yield read info
        yield isize, r1.is_reverse
        if (max_frags is not None) and (frags >= max_frags):
            break
    logging.debug("Processed fragments: %d" % (frags))
    logging.debug("Unmapped fragments: %d" % (unmapped))

def inspect_random(mean, stdev, min_isize, max_isize, max_frags=100000):
    frags = 0
    unmapped = 0
    while True:
        frags += 1
        isize = int(round(random.normalvariate(mean, stdev),0))
        if not (min_isize <= isize <= max_isize):
            unmapped += 1
            continue
        is_reverse = bool(random.randint(0,1))
        yield isize, is_reverse
    logging.debug("Processed fragments: %d" % (frags))
    logging.debug("Unmapped fragments: %d" % (unmapped))

class RnaseqLibraryCharacteristics(object):
    def __init__(self, min_frag_size, max_frag_size):
        self.min_frag_size = max(0, min_frag_size)
        self.max_frag_size = max(0, max_frag_size)
        self.arr = array.array('L', (0 for x in xrange(self.min_frag_size, self.max_frag_size+1)))
        self.read1_rev_count = 0
        self.read1_count = 0

    def from_stream(self, myiter):
        total = 0
        inrange = 0
        for isize, is_reverse in myiter:
            if isize is not None:
                if (self.min_frag_size <= isize <= self.max_frag_size):
                    self.arr[isize - self.min_frag_size] += 1
                    inrange += 1
            self.read1_rev_count += int(is_reverse)
            self.read1_count += 1
            total += 1
        logging.debug("Stream had %d total samples with %d within the fragment size range" % (total, inrange))

    def strandedness(self):
        return self.read1_rev_count / float(self.read1_count)

    def isize_at_percentile(self, per):
        n = sum(self.arr)
        per_n = n * per / 100.0
        count = 0
        for isize,x in enumerate(self.arr): 
            count += x
            if (count >= per_n):
                break
        return isize + self.min_frag_size
    
    def percentile_at_isize(self, isize):
        if isize < self.min_frag_size:
            return 0.0
        elif isize > self.max_frag_size:
            return 100.0
        ind = isize - self.min_frag_size
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
            isize = i + self.min_frag_size
            std = std + x*((isize - mean)**2)
            n += x
        std = (std / float(n-1))**0.5
        return std

    def plot_frag_size_dist(self, filename, max_percentile=100.0):
        import matplotlib.pyplot as plt
        plot_max_isize = self.isize_at_percentile(max_percentile)
        xdata = range(self.min_isize, plot_max_isize+1)
        num_samples = self.num_frag_size_samples()
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
        print >>fileh, '\t'.join(["#read1_rev_count", "read1_count"])
        print >>fileh, '\t'.join(map(str, [self.read1_rev_count, self.read1_count]))
        print >>fileh, '\t'.join(["#insert_size", "num_samples"])
        for i,x in enumerate(self.arr):
            print >>fileh, '\t'.join([str(i + self.min_isize), str(x)])        

    @staticmethod
    def from_file(fileh):
        res = RnaseqLibraryCharacteristics()
        # skip first comment
        fileh.next()
        # read strandedness
        fields = map(int, fileh.next().strip().split('\t'))
        res.read1_rev_count = fields[0]
        res.read1_count = fields[1]
        # read frag size distribution
        isizes = []
        counts = []
        for line in fileh:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            i,x = map(int, fields[0:2])
            isizes.append(i)
            counts.append(x)
        res.min_isize = isizes[0]
        res.max_isize = isizes[-1]
        res.arr = array.array('L', counts)
        return res 


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

