'''
Contains classes for generating text output for a Sample object.

@author: cnielsen
'''
import logging, datetime

logger = logging.getLogger(__name__)

class TextWriter:
    
    def generateTextOutput(self, version, sample, ofilename, ):
        istring = "writing text output to file " + ofilename
        logger.info(istring)
        print istring
        
        fh = open(ofilename, 'w')
        
        # header
        fh.write("#museqreportbuilder.py v" + version + " " + str(datetime.datetime.now()) + "\n")
        
        fh.write("\n#probability density")
        pDensity = sample.getProbDensity()
        for fType in sorted(pDensity.keys()):
            fh.write("\n#" + fType + "\n")
            fh.write("prob\tcount\n")
            for p in sorted(pDensity[fType].keys()):
                fh.write(str(p) + "\t" + str(pDensity[fType][p]) + "\n")
                
        thresholds = sample.getThresholds()
        fh.write("\n#probability versus allele ratio (P > " + str(thresholds[-1]) + ")\n")
        fh.write("prob\tallele ratio\tcount\n")     
        (probs, ratios, counts) = sample.getAlleleRatioAndProbWithCounts(thresholds[-1])
        for i,p in enumerate(probs):
            fh.write(str(p) + "\t" + str(ratios[i]) + "\t" + str(counts[i]) + "\n")
        
        for threshold in thresholds:
            fh.write("\n#coverage versus allele ratio (P > " + str(threshold) + ")\n")  
            fh.write("cov\tallele ratio\tcount\n")     
            (covs, ratios, counts) = sample.getAlleleRatioAndCovWithCounts(threshold)
            for i,c in enumerate(covs):
                fh.write(str(c) + "\t" + str(ratios[i]) + "\t" + str(counts[i]) + "\n")
                        
        if sample.hasTrinucData():
            fh.write("\n#trinucleotide contexts (P > " + str(thresholds[-1]) + ")\n")
            fh.write("mut\tcontext\tcount\n")
            (mutations, contextSets, countSets) = sample.getTriNucData()
            for i,m in enumerate(mutations):
                for j,c in enumerate(contextSets[i]):
                    fh.write(m + "\t" + c + "\t" + str(countSets[i][j]) + "\n")
                
        fh.write("\n#genomic distribution of mutations (P > " + str(thresholds[-1]) + ")\n")
        fh.write("chrom\tbin start\tmut per Mbp\n")
        posDensitySets = sample.getGlobalPosDensity()
        cSizes = sample.getChromSizes()
        for i,c in enumerate(sample.getChromNames()):
            bSize = cSizes[i]/len(posDensitySets[i])
            for j,d in enumerate(posDensitySets[i]):
                fh.write(c + "\t" + str(int(j * bSize)) + "\t" + str(d) + "\n")
        fh.close()
        logger.info("text output complete")
        