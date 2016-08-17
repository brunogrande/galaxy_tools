'''
Parser for mutationSeq VCF output

@author: cnielsen
'''
import logging
from os.path import basename

import museqsample

logger = logging.getLogger(__name__)

class VCFParser:
    
    def processVCFInfo(self, iString):
        labels = ["PR", "TR", "TA", "NR", "NA", "TC"]
        # tolerate these variations in field names
        possibleAltLabels = ["TUMOUR_REF", "TUMOUR_ALT", "NORMAL_REF", "NORMAL_ALT"] 
        values = [None] * len(labels)
        # parse the INFO string key-value pairs
        for field in iString.split(";"):
            k, v = [f.strip() for f in field.split("=")]
            if k in labels:
                # all values are floats except for the trinucleotide context "TC"
                if v == ".":
                    v = None
                elif k != "TC" and v != ".":
                    v = float(v)
                values[labels.index(k)] = v
            else:
                # handle possible alternatives
                if k in possibleAltLabels:
                    # e.g. "NORMAL_REF=C,39"
                    v = float(v.split(",")[1])
                    if k == "TUMOUR_REF":
                        values[labels.index("TR")] = v
                    elif k == "TUMOUR_ALT":
                        values[labels.index("TA")] = v
                    elif k == "NORMAL_REF":
                        values[labels.index("NR")] = v
                    elif k == "NORMAL_ALT":
                        values[labels.index("NA")] = v
        return values
    
    def processVCFLine(self, iString):
        fields = iString.strip().split("\t")
        baseValues = fields[:7]
        
        # adjust types
        pIndex = 1; qIndex = 5
        for i,f in enumerate(baseValues):
            if f == ".":
                baseValues[i] = None
            elif i == pIndex:
                baseValues[pIndex] = int(baseValues[pIndex]) # positions
            elif i == qIndex:
                baseValues[qIndex] = float(baseValues[qIndex]) # quality values

        # process INFO string
        infoValues = self.processVCFInfo(fields[7])
        
        return baseValues + infoValues
            
    def parseSample(self, fn, thresholds, fType, build):
        """
        Parses input VCF file into a Sample object with data from input file
        """
        istring = "parsing " + fn
        logger.info(istring)   
        print istring
        
        sample = museqsample.Sample(thresholds, fType, build, basename(fn))
        fh = open(fn, 'r')
        count = 0
        for l in fh.readlines():
            if l.startswith("#"):
                # logger.info("skipping header line - " + l.strip())
                continue
            try:
                (c, pos, id, r, a, q, f, pr, tr, ta, nr, na, tc) = self.processVCFLine(l)
                sample.addValue(filter_type = fType, chrom = c, pos = pos, ref = r, alt = a, qual = q, filter = f, 
                                prob=pr, ref_count = tr, alt_count = ta, 
                                norm_ref_count = nr, norm_alt_count = na, trinuc = tc)
                count += 1
            except Exception,e:
                logger.error(str(e) + "\nline: '" + l + "'")
                raise e
        fh.close()
        logger.info("skipped cases with unknown nucleotide: " + str(sample.getUnknownNucCount()))
        logger.info("parsing complete (" + str(count) + " mutations)")
        logger.info("filter counts: " + sample.getStats())

        return sample