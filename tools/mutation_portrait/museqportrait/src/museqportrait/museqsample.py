'''
Contains classes for managing mutationSeq output data

@author: cnielsen
'''
import logging
import random
from collections import defaultdict

logger = logging.getLogger(__name__)

class NucleotideError(Exception):
    def __init__(self, seq):
        self.seq = "Unknown nucleotide sequence: " + seq
    def __str__(self):
        return self.seq

class Sample:
       
    CHROM_NAMES = [str(i+1) for i in range(22)] + ['X', 'Y']
    HG18_CHROM_SIZES = [247249719.0, 242951149.0, 199501827.0, 191273063.0, 180857866.0, 170899992.0, 158821424.0, 146274826.0,
                        140273252.0, 135374737.0, 134452384.0, 132349534.0, 114142980.0, 106368585.0, 100338915.0, 88827254.0,
                        78774742.0, 76117153.0, 63811651.0, 62435964.0, 46944323.0, 49691432.0, 154913754.0, 57772954.0]
    
    HG19_CHROM_SIZES = [249250621.0, 243199373.0, 198022430.0, 191154276.0, 180915260.0, 171115067.0, 159138663.0, 146364022.0,
                        141213431.0, 135534747.0, 135006516.0, 133851895.0, 115169878.0, 107349540.0, 102531392.0, 90354753.0,
                        81195210.0, 78077248.0, 59128983.0, 63025520.0, 48129895.0, 51304566.0, 155270560.0, 59373566.0]
    
    def __init__(self, thresholds, fType, build, name):
        if build == "hg19":
            self.CHROM_SIZES = self.HG19_CHROM_SIZES
        elif build == "hg18":
            self.CHROM_SIZES = self.HG18_CHROM_SIZES
        else:
            eString = "Unknown genome build: " + build
            e = Exception(eString)
            logger.error(eString)
            raise e
        self.NUM_CHROM_BINS = 3000 # results in overplotting, but better for resolution of frequency spikes; produces bins of ~1 Mbp
        self.chromBinCounts, self.chromBinSizes = self.computeChromBinSize(self.NUM_CHROM_BINS) # assign number of bins in proportion to chrom size
        
        # labels and stats
        self.name = name
        self.filterType = fType
        self.thresholds = thresholds
        self.stats = defaultdict(lambda : defaultdict(int))
        # self.filterStats = defaultdict(int) # mutation counts by different filter types
        # self.thresholdStats = defaultdict(int) # mutation counts by different probability thresholds
        self.unknown_nuc_count = 0 # reported in log
        
        # mutation counts across a range of probabilities separated by FILTER type (ALL, PASS, INDEL, etc)
        self.probDensity = defaultdict(lambda : defaultdict(int))
        
        # dictionary for each threshold - each storing [coverage][allele_ratio] = counts
        self.coverageToAlleleRatio = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))
        # dictionary for each threshold - each storing [probability][allele_ratio] = counts
        self.probToAlleleRatio = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))

        # mutation counts by tri-nucleotide types
        self.triNucDensity = defaultdict(lambda : defaultdict(int))
        
        # mutation counts across genomic increments
        self.posDensity = self.initializeChromosomes()
        # list of lists of normalized mutation counts per chromosome 
        self.globalPosDensity = None

        # titan results (list of tuples) indexed by chromosome: dict[chrom] = [(position, logRatio, allelicRatio, copyNumber, copyNumberCall), ...]
        self.TitanResults = defaultdict(list)
        self.POSITION_INDEX = 0
        self.LOG_RATIO_INDEX = 1
        self.ALLELIC_RATIO_INDEX = 2
        self.COPY_NUMBER_INDEX = 3
        self.COPY_NUMBER_CALL_INDEX = 4
        
        # genome-wide titan results (list of tuples) sampled by different variables 
        self.sampledByPosAndLogRatio = None # for storing subsampled data
        self.sampledByPosAndAllelicRatio = None # for storing subsampled data
        
    def addTriNuc(self, mut, tri_context):
        try:
            self.triNucDensity[mut][tri_context] += 1
        except:
            e = NucleotideError(mut + " or " + tri_context)
            logger.warning(e)
            raise e
    
    def addValue(self, **options):
        try:
            # filter based on lowest threshold
            t_lower = self.thresholds[0]
            if options['prob'] < t_lower: return
            if options['filter_type'] != None and options['filter'] != options['filter_type']: return
            
            self.stats[t_lower]["ALL"] += 1
            if options['filter'] != None: 
                self.probDensity[options["filter"]][options["prob"]] += 1
                self.stats[t_lower][options["filter"]] += 1
                
            self.probDensity["ALL"][options["prob"]] += 1
            
            coverage = options["alt_count"] + options["ref_count"]
            allele_ratio = options["alt_count"]/coverage
            self.coverageToAlleleRatio[t_lower][coverage][allele_ratio] += 1
            self.probToAlleleRatio[t_lower][options["prob"]][allele_ratio] += 1
            
            # record counts and other stats above other specified thresholds
            for t in self.thresholds[1:]:
                if options['prob'] >= t:
                    self.stats[t]["ALL"] += 1
                    if options['filter'] != None: 
                        self.stats[t][options['filter']] += 1
                    self.coverageToAlleleRatio[t][coverage][allele_ratio] += 1
                    self.probToAlleleRatio[t][options["prob"]][allele_ratio] += 1
                
            # store trinuc and genomic position for only highest threshold values    
            if options["prob"] >= self.thresholds[-1]:
                if options['trinuc'] != None:
                    self.storeTriNucContext(options['ref'], options['alt'], options['trinuc'])
                self.posDensity[options['chrom']][options['pos']] += 1
        except Exception,e:
            logger.error(e)
            raise e
            
    def addCopyNumber(self, chrom, position, logRatio, allelicRatio, copyNumber, callType):
        self.TitanResults[chrom].append([position, logRatio, allelicRatio, copyNumber, callType])

    def assembleValues(self, dictOfDict):
        # can have multiple dependent values for each independent value
        iValues = []; dValues = []; 
        for i in sorted(dictOfDict.keys()):
            for d in dictOfDict[i].keys():
                count = dictOfDict[i][d]
                dValues.extend([d]*count)
                iValues.extend([i]*count)
        return (iValues, dValues)
    
    def assembleValuesWithCounts(self, dictOfDict):
        # can have multiple dependent values for each independent value
        iValues = []; dValues = []; counts = []
        for i in sorted(dictOfDict.keys()):
            for d in dictOfDict[i].keys():
                dValues.append(d)
                iValues.append(i)
                counts.append(dictOfDict[i][d])
        return (iValues, dValues, counts)   
    
    def computeChromBinSize(self, n):
        # assign number of bins in proportion to chrom size
        genomeLen = sum(self.CHROM_SIZES)
        bCounts = []; bSizes = [] 
        for size in self.CHROM_SIZES:
            f = float(size)/genomeLen
            bCount = int(round(f * n))
            bCounts.append(bCount)
            bSizes.append(size/bCount)
        return (bCounts, bSizes)

    def computeGlobalPosDensity(self):
        self.globalPosDensity = []
        # sum counts for each chrom - convert to mutations/Mbp
        for i,chrom in enumerate(self.CHROM_NAMES):
            positions = sorted(self.posDensity[chrom].keys())
            counts = [0]*self.chromBinCounts[i]; 
            j = 0; threshold = self.chromBinSizes[i]
            for pos in positions:
                while pos > threshold and j < len(counts)-1:
                    j += 1
                    threshold += self.chromBinSizes[i]
                counts[j] += 1
            # normalize by bin size and store as counts per Mbp
            counts = [(c*10**6)/self.chromBinSizes[i] for c in counts]
            self.globalPosDensity.append(counts)
    
    def getAlleleRatioAndCov(self, threshold):
        return self.assembleValues(self.coverageToAlleleRatio[threshold])
    
    def getAlleleRatioAndCovWithCounts(self, threshold):
        return self.assembleValuesWithCounts(self.coverageToAlleleRatio[threshold])
    
    def getAlleleRatioAndProb(self, threshold):
        return self.assembleValues(self.probToAlleleRatio[threshold])
    
    def getAlleleRatioAndProbWithCounts(self, threshold):
        return self.assembleValuesWithCounts(self.probToAlleleRatio[threshold])
    
    def getChromBinCounts(self):
        return self.chromBinCounts
    
    def getChromBinSize(self, chrom):
        s = None
        for i,c in enumerate(self.CHROM_NAMES):
            if c == chrom:
                s = self.chromBinSizes[i]
                break
        return s
    
    def getChromNames(self):
        return self.CHROM_NAMES
    
    def getChromSizes(self):
        return self.CHROM_SIZES
    
    def getFilterType(self):
        return self.filterType
    
    def randomlySample(self, iTupleSet):
        if len(iTupleSet) > 1: 
            i = random.sample(xrange(len(iTupleSet)), 1)[0]
            return iTupleSet[i]
        elif len(iTupleSet) == 1:
            return iTupleSet[0]
        else:
            eString = "Cannot randomly sample from an empty set"
            e = Exception(eString)
            logger.error(e)
            raise e
    
    def sampleTupleSetByVar(self, iTupleSet, tIndex, varMin, varBinSize):
        sampledSet = []
        if len(iTupleSet) > 1:
            # sort the tuples by var at index tIndex
            sorted_by_var = sorted(iTupleSet, key=lambda tup: tup[tIndex])
            varThres = varMin + varBinSize
            currTupleSet = []
            for tup in sorted_by_var:
                while tup[tIndex] > varThres:
                    # process current set if any
                    if len(currTupleSet) > 0:
                        sampledSet.append(self.randomlySample(currTupleSet))
                    # reset
                    currTupleSet = []
                    varThres += varBinSize
                # store current tuple for later processing
                currTupleSet.append(tup)
            if len(currTupleSet) > 0:
                # process last bin if any
                sampledSet.append(self.randomlySample(currTupleSet))
        elif len(iTupleSet) == 1:
            # nothing to sample
            sampledSet = iTupleSet
        return sampledSet
    
    def sampleTupleSetByTwoVars(self, iTupleSet, tIndex1, tIndex2, varMin1, varMin2, varBinSize1, varBinSize2):
        sampledSet = []
        if len(iTupleSet) > 1:
            # sort the tuples by first var at index tIndex1
            sorted_by_var_1 = sorted(iTupleSet, key=lambda tup: tup[tIndex1])
            varThres1 = varMin1 + varBinSize1
            currTupleSet = []
            for tup in sorted_by_var_1:
                while tup[tIndex1] > varThres1:
                    # process current set if any
                    if len(currTupleSet) > 0:
                        # sample current tuple set according to second variable at index tIndex2
                        sampledSet.extend(self.sampleTupleSetByVar(currTupleSet, tIndex2, varMin2, varBinSize2))
                    # reset
                    currTupleSet = []
                    varThres1 += varBinSize1
                # store current tuple for later processing
                currTupleSet.append(tup)
            if len(currTupleSet) > 0:
                # process last bin if any
                sampledSet.extend(self.sampleTupleSetByVar(currTupleSet, tIndex2, varMin2, varBinSize2))
        elif len(iTupleSet) == 1:
            # nothing to sample
            sampledSet = iTupleSet
        return sampledSet       
        
    def sampleTitanGenomeWide(self, gBinCount, varRange, varBinCount):
        # bin and sample genome-wide and across the specified var
        sampledTuples = []   
        binCounts, binSizes = self.computeChromBinSize(gBinCount)
        for i,chrom in enumerate(self.CHROM_NAMES):
            currSampledTuples = self.sampleTitanChromWide(chrom, self.CHROM_SIZES[i], binCounts[i], varRange, varBinCount)
            if len(currSampledTuples) > 0:
                # update genomic positions from chrom-centric to global
                currSampledTuplesAdjusted = [] # ensures a copy
                baseSize = sum(self.CHROM_SIZES[:i])
                for i in range(len(currSampledTuples)):
                    currTuple = list(currSampledTuples[i]) # copy
                    currTuple[self.POSITION_INDEX] += baseSize
                    currSampledTuplesAdjusted.append(currTuple)
                    if currSampledTuplesAdjusted[i][self.POSITION_INDEX] > 491986559.0:
                        print "error: after update longer than threshold " + str(currSampledTuplesAdjusted[i])
                sampledTuples.extend(currSampledTuplesAdjusted)
        return sampledTuples
            
    def sampleTitanChromWide(self, chrom, chromSize, gBinCount, varRange, varBinCount):
        sampledTuples = []
        if self.TitanResults.has_key(chrom):
            # compute var bin sizes
            gBinSize = int(chromSize/gBinCount)
            varMin = varRange[0]
            varBinSize = (varRange[1] - varRange[0]) / varBinCount
            # sample and bin
            sampledTuples = self.sampleTupleSetByTwoVars(self.TitanResults[chrom], self.POSITION_INDEX, self.LOG_RATIO_INDEX, 0, varMin, gBinSize, varBinSize)
        return sampledTuples
    
    def getAllelicRatiosGenomeWide(self, subsample, gBinCount, arRange, arBinCount, callType):
        if subsample:
            if self.sampledByPosAndAllelicRatio == None:
                print "loh " + str(gBinCount) + " " + str(arRange) + " " + str(arBinCount)
                self.sampledByPosAndAllelicRatio = self.sampleTitanGenomeWide(gBinCount, arRange, arBinCount)
            subset = [tup for tup in self.sampledByPosAndAllelicRatio if tup[self.COPY_NUMBER_CALL_INDEX] == callType]
            xPositions = [tup[self.POSITION_INDEX] for tup in subset]
            logRatios = [tup[self.ALLELIC_RATIO_INDEX] for tup in subset]
            return xPositions, logRatios
        else:
            # TODO - return all
            return None      
    
    def getLogRatiosGenomeWide(self, plot_all, gBinCount, lrRange, lrBinCount, copyType):
        if plot_all:
            # TODO - return all
            return None
        else:
            if self.sampledByPosAndLogRatio == None:
                print "copy number " + str(gBinCount) + " " + str(lrRange) + " " + str(lrBinCount)
                self.sampledByPosAndLogRatio = self.sampleTitanGenomeWide(gBinCount, lrRange, lrBinCount)
            subset = [tup for tup in self.sampledByPosAndLogRatio if tup[self.COPY_NUMBER_INDEX] == copyType]
            xPositions = [tup[self.POSITION_INDEX] for tup in subset]
            logRatios = [tup[self.LOG_RATIO_INDEX] for tup in subset]
            return xPositions, logRatios
    
    def getName(self):
        return self.name

    def getGlobalPosDensity(self):
        if self.globalPosDensity == None:
            self.computeGlobalPosDensity()
        return self.globalPosDensity

    def getProbDensity(self):
        return self.probDensity
    
    def getTriNucData(self):
        mutations, contextSets, countSets = ([], [], [])
        # construct all possible contexts - some may not have been observed in self.triNucDensity
        pyr = ['C', 'T']
        base = ['A', 'C', 'G', 'T']
        for p in pyr:
            for b in base:
                if p != b:
                    mut = p + ">" + b
                    mutations.append(mut)
                    tri_contexts = [base[i] + p + base[j] for i in range(4) for j in range(4)]
                    counts = [0]*len(tri_contexts)
                    if self.triNucDensity.has_key(mut):
                        for i,c in enumerate(tri_contexts):
                            if self.triNucDensity[mut].has_key(c):
                                counts[i] = self.triNucDensity[mut][c]
                    contextSets.append(tri_contexts)
                    countSets.append(counts)
        return (mutations, contextSets, countSets)
    
    def getThresholds(self):
        return self.thresholds
    
    def getStats(self):
        s = ""
        for thres in self.stats.keys():
            s += 'prob > ' + str(thres) + ": " + str(self.stats[thres]["ALL"])
            if len(self.stats[thres].keys()) > 1:
                s += " ("
                s += "; ".join(['{0}: {1}'.format(t,c) for (t,c) in self.stats[thres].items() if t != "ALL"])
                s += ")"
            s += "\n"
        return s
    
    def getUnknownNucCount(self):
        "Number of nucleotides encountered that are not 'A', 'C', 'G', or 'T'"
        return self.unknown_nuc_count
    
    def hasTrinucData(self):
        if self.triNucDensity.keys() != []:
            return True
        else:
            return False
    
    def hasCopyNumberData(self):
        if self.TitanResults.keys() != []:
            return True
        else:
            return False
    
    def initializeChromosomes(self):
        d = {}
        for s in self.CHROM_NAMES:
            d[s] = defaultdict(int)
        return d
        
    def reverseComplement(self, s):
        c = ""
        for letter in s:
            if letter == "A":
                c += "T"
            elif letter == "T":
                c += "A"
            elif letter == "C":
                c += "G"
            elif letter == "G":
                c += "C"
            else:
                e = NucleotideError(s)
                logger.warning(e)
                raise e
        c = c[::-1] # reverse the string
        return c
                
    def storeTriNucContext(self, ref, alt, tri_ref):
        try:
            if (ref == 'A' or ref == 'G'):
                ref = self.reverseComplement(ref)
                alt = self.reverseComplement(alt)
                tri_ref = self.reverseComplement(tri_ref)                   
            self.addTriNuc(ref + ">" + alt, tri_ref)
        except NucleotideError:
            # these errors have been logged
            self.unknown_nuc_count += 1