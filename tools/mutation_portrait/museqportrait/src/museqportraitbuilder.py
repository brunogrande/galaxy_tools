"""
Script for generating a summary report of mutationSeq output

@author: cnielsen
"""
import argparse
import logging
import os
import traceback

from museqportrait import VCFParser, TitanParser, Plotter, TextWriter

VERSION = "0.99.0"
logger = logging.getLogger(__name__)

def parse_args():
    desc = "This script generates a 'portrait' summarizing mutationSeq output. " + \
           "Input can either be one or more comma separated VCF files or a " + \
           "directory name. If a directory name is provided, all VCF files " + \
           "contained in that directory will be processed. One 'portrait' pdf/" + \
           "text file pair is generated for each VCF file."
    argparser = argparse.ArgumentParser(description=desc)
    # required arguments
    argparser.add_argument('input', metavar='INPUT', nargs='+',
                           help='one or more VCF files or directory containing VCF files')    
    # optional arguments
    lowerDefault = '0.5'
    upperDefault = '0.9'
    buildDefault = 'hg19'    
    argparser.add_argument('--lower', default=lowerDefault, 
                           help='lower probability threshold (default: ' + lowerDefault + ")")

    argparser.add_argument('--upper', default=upperDefault, 
                           help='upper probability threshold (default: ' + upperDefault + ")")
    argparser.add_argument('--filter', default=None, 
                           help='include only those mutations with this filter state (e.g. PASS, INDEL, etc.)')
    argparser.add_argument('--build', default=buildDefault,
                           help='genome build to use for genome-wide plots (default: ' + buildDefault + '; ' + \
                           'only hg18 and hg19 currently supported)')
    # Copy number code needs further debugging - reserving for next release
    # argparser.add_argument('--cnv', default=None, 
    #                        help='name of TITAN copy number results file')
    # argparser.add_argument('--ploidy', default=None,
    #                        help='estimated ploidy reported by TITAN')
    # argparser.add_argument('--plot-all', action='store_true',
    #                        help='There are often more data points than can be visually distinguished in a plot. ' + \
    #                        'By default, points will be sub-sampled. Use this flag to override the sub-sampling and ' + \
    #                        'plot all points (WARNING: may result in large output files)')
    argparser.add_argument('--out-dir', default=".", 
                           help='name of output directory in which to write output files (.pdf and .txt)')
    argparser.add_argument('--version', action="version", version=VERSION)
    return argparser.parse_args()
       
def process_input(inputList):
    files = []; dirs = []
    for s in inputList:
        if os.path.isfile(s):
            if is_vcf(s):
                files.append(s)
        elif os.path.isdir(s):
            dirs.append(s)
        else:
            eString = "Input is neither a file nor a directory: " + s
            e = Exception(eString)
            logger.error(eString)
            raise e
        
    # get all files in specified dirs (if any)
    for d in dirs:
        for f in os.listdir(d):
            ifn = os.path.join(d, f)
            if os.path.isfile(ifn) and is_vcf(ifn):
                files.append(ifn)
    return files

def process_output(ifiles, odir, thres, ftype):
    oplotfiles = []; otextfiles = []
    if os.path.isdir(odir):
        for p in ifiles:
            ifn = ".".join(os.path.basename(p).split(".")[:-1])
            ofn = "_".join([ifn, str(thres[0]), str(thres[-1])])
            if ftype != None: 
                ofn = ofn + "_" + ftype
            oplotfile = ofn + ".pdf"
            otextfile = ofn + ".txt"
            oplotfiles.append(os.path.join(odir, oplotfile))
            otextfiles.append(os.path.join(odir, otextfile))
    else:
        estring = "Value of outdir must be a directory: " + odir
        e = Exception(estring)
        logger.error(estring)
        raise e
    return oplotfiles, otextfiles

def is_vcf(fName):
    if fName.split(".")[-1].lower() == "vcf":
        return True
    else:
        return False    
            
def main():
    try:
        logging.basicConfig(filename="museqportraitbuilder.log", filemode="w", level=logging.DEBUG, 
                            format='%(asctime)s - %(levelname)s - %(name)s - %(funcName)s - line %(lineno)d - %(message)s')
        logger.info("begin")
        
        # collect user inputs and defaults
        args = parse_args()
        thresholds = (float(args.lower), float(args.upper))
        # handles cases where input name is a directory
        ifilenames = process_input(args.input)
        # auto-generate output file names
        oplotfiles, otextfiles = process_output(ifilenames, args.out_dir, thresholds, args.filter)
        
        # temporary while improving copy-number plotting code
        args.cnv = None
        args.plot_all = False
        
        # generate one report per input VCF file
        vcfparser = VCFParser(); 
        titanparser = TitanParser();
        plotter = Plotter(); 
        writer = TextWriter();
        sample = None
        for i,fn in enumerate(ifilenames):
            sample = vcfparser.parseSample(fn, thresholds, args.filter, args.build)
            if (args.cnv):
                titanparser.parseTitanResultsFile(args.cnv, sample, args.ploidy)
            plotter.generateSampleReport(sample, oplotfiles[i], args.plot_all)
            writer.generateTextOutput(VERSION, sample, otextfiles[i])
        logger.info("complete")
    except Exception, e:
        logger.error(e)
        print e
        print traceback.print_exc()
    
if __name__ == "__main__":
    main()