'''
Parser for Titan copy number output

@author: cnielsen
'''
import logging
import math

logger = logging.getLogger(__name__)

class TitanParser:
    
    def __init__(self):
        # ignore segments that are shorter than this threshold
        self.MIN_SEG_LENGTH = 5000
    
    def convertCNtoType(self, cn):
        """No longer being used"""
        cType = ""
        if (cn == 0):
            cType = "HOMD"
        elif (cn == 1):
            cType = "DLOH"
        elif (cn == 2):
            cType = "NEUT"
        elif (cn == 3):
            cType = "GAIN"
        elif (cn == 4 or cn == 5):
            cType = "AMP"
        else:
            eString = "Unexpected copy number: " + cn
            e = Exception(eString)
            logger.error(eString)
            raise e
        return cType       
    
    def parseTitanResultsFile(self, fn, sample, ploidy):
        istring = "parsing " + fn
        logger.info(istring)   
        print istring
        fh = open(fn, 'r')
        
        for l in fh.readlines():
            if l.startswith("Chr"): continue
            try:
                fields = l.strip().split("\t")
                chrom = fields[0]
                position = float(fields[1]) # use float to handle case of scientific notation
                alleleRatio, logRatio = [float(v) for v in fields[5:7]]
                if ploidy != None: logRatio += math.log(float(ploidy)/2, 2) # correct for ploidy
                copyNumber = fields[7]
                titanCall = fields[9]
                sample.addCopyNumber(chrom, position, logRatio, alleleRatio, copyNumber, titanCall)                
            except Exception,e:
                logger.error(str(e) + "\nline: '" + l + "'")
                raise e
        fh.close()
        logger.info("copy number parsing complete")
    
    def parseTitanSegFile(self, fn, sample):
        """No longer being used"""
        istring = "parsing " + fn
        logger.info(istring)   
        print istring
        fh = open(fn, 'r')
        for l in fh.readlines():
            if l.startswith("Sample"): continue
            try:
                fields = l.strip().split("\t")
                l = int(fields[4])
                if l < self.MIN_SEG_LENGTH: continue
                center = int(fields[2]) + l/2
                cType = self.convertCNtoType(int(fields[9]))
                sample.addCopyNumberValue(fields[1], center, float(fields[6]), cType)
            except Exception,e:
                logger.error(str(e) + "\nline: '" + l + "'")
                raise e
        fh.close()
        logger.info("copy number parsing complete")