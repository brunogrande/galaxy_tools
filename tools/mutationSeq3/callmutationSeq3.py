#!/usr/bin/env python

import sys
import os
import getopt
import subprocess
import fileHandler


#python classify.py normal:/data3/rdmorin/public/dmestry_data/A01458_9_lanes_dupsFlagged.bam tumour:/data3/rdmorin/public/dmestry_data/A01456_10_lanes_dupsFlagged.bam reference:/home/dmestry/data/human_all.fasta model:model.npz --config metadata.config --interval 9 --threshold 0.001 --out /home/dmestry/outputs/0725Chr9mutationSeqOut1.vcf


# run readCounter on input1
def runClassify(normalBam, tumourBam,refFasta,model,config,output,interval=None,threshold=None):
	
	#change classifyDir to where mutaitonSeq 3.0.0 is installed.
	classifyDir = "/home/bgrande/software/mutationSeq_3.0.1"
	
	classifyCommand = ["python ", classifyDir + "/classify.py"," normal:" + normalBam, " tumour:" + tumourBam, " reference:" + refFasta, " model:" + model, " --config" + ' ' + config, " --out " + output]
	
	if interval:
		classifyCommand.append(" --interval " + interval)
		print "added interval option"
	else:
		print "no interval"
	
	if threshold:
		classifyCommand.append(" --threshold " + threshold)
		print "added threshold option"
	else:
		print "no threshold"
	
	print classifyCommand
	print
		
	cString=""
	
	joinedcString = cString.join(classifyCommand)
	print joinedcString
	subprocess.check_call(joinedcString, shell=True, stderr=subprocess.STDOUT)


# copy galaxy bam index from galaxy directory to data directory
def copyGalaxyBamIndexToDataDir(galaxyDir,dataDir,galaxyBamIndexFile,bamFileName):
	copyCommand=['cp', galaxyDir + '/' + galaxyBamIndexFile, dataDir + '/' + bamFileName + '.bai']
	print "copy command is "
	print copyCommand
	subprocess.call(copyCommand)		


#take file, including path, return path and file 
def splitArgs(arg):
	path, File = fileHandler.getFileAndPath(arg)
	return path, File


def main():
	#get options and args
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hn:t:o:r:c:i:m:", ["threshold=","interval=","normalBamIndex=","tumourBamIndex="])
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		usage()
		sys.exit(2)
	
	output = None
	normalBam = None
	tumourBam = None
	refFasta = None
	config = None
	interval = None
	model = None
	threshold = None
	tumourBamIndex = None
	normalBamIndex = None
	
	#parse options and arguements
	for option, arg in opts:
		if option == "-v":
			continue
		elif option in ("-h", "--help"):
			usage()
			sys.exit()
		elif option in ("-o", "--output", "--out"):
			output = arg
		elif option in ("-n", "--normal"):
			normalBam = arg
		elif option in ("-t", "--tumour"):
			tumourBam = arg
		elif option in ("--tumourBamIndex"):
			tumourBamIndex = arg
		elif option in ("--normalBamIndex"):
			normalBamIndex = arg
		elif option in ("-r", "--reference"):
			refFasta = arg
		elif option in ("-c", "--config"):
			config = arg
		elif option in ("-i", "--interval"):
			interval = arg
		elif option in ("-m", "--model"):
			model = arg
		elif option in ("--threshold"):
			threshold = arg
		else:
			assert False, "unhandled option"
	
	
	#what the model and config files are
	config = "/home/bgrande/software/mutationSeq_3.0.1/metadata.config"
	model = "/home/bgrande/software/mutationSeq_3.0.1/model.npz"
	
	
	#print all inputs
	normal_data_dir,normal_bamFile_name = splitArgs(normalBam)
	tumour_data_dir,tumour_bamFile_name = splitArgs(tumourBam)
	
	print normal_data_dir,normal_bamFile_name 
	print tumour_data_dir,tumour_bamFile_name 
	print
	galaxy_dir,galaxy_bamIndex = splitArgs(normalBamIndex)
	galaxy_dir2,galaxy_bamIndex2 = splitArgs(tumourBamIndex)
	
	print galaxy_dir,galaxy_bamIndex
	print galaxy_dir2,galaxy_bamIndex2
	print
	print "reference fasta is, ", refFasta
	
	print "config is, ", config
	print "interval is, ", interval
	print "model is, ", model
	print "threshold is, ", threshold
	print
	print "output file is, ", output
	print
	
	
	#execute supbrocesses to run classify
	
	#for normal bam file
	copyGalaxyBamIndexToDataDir(galaxy_dir,normal_data_dir,galaxy_bamIndex,normal_bamFile_name)
	
	#for tumour bam file
	copyGalaxyBamIndexToDataDir(galaxy_dir2,tumour_data_dir,galaxy_bamIndex2,tumour_bamFile_name)
	
	#run mutationSeq classify
	if interval and threshold:
		runClassify(normalBam,tumourBam,refFasta,model,config,output,interval=interval,threshold=threshold)
	elif threshold:
		runClassify(normalBam,tumourBam,refFasta,model,config,output,threshold=threshold)
	elif interval:
		runClassify(normalBam,tumourBam,refFasta,model,config,output,interval=interval)
	else:
		runClassify(normalBam,tumourBam,refFasta,model,config,output)
	
	print "done"




if __name__ == "__main__":
	main()

