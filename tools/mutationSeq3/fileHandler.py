import sys
import os.path


def openFile(pathFile):
	try:
		In = open(pathFile,"r+")
		return In 
	except IOError:
		print "couldn't read input file "


def createOutputFile(output):
	try:
		out = open(output,"a+")
		return out
	except IOError:
		print "couldn't create output file"

def processMultipleFiles(directory):
	files = []
	#get file names from directory
	for dirname, dirnames, filenames in os.walk(directory):
		# print path to all filenames
	    for filename in filenames:
			#if filename.endswith(".txt"):
			files.append(os.path.join(dirname, filename))
	return files	


def enterPath():
	#asks user to enter file with path
	try:
		f = raw_input("Please enter a path:")
		return f
	except:
		print "include the path"


def getFileName(f):
	#gets filename from the path entered
	fname = os.path.basename(f)
	return fname	


def getFileAndPath(fp):
	path, File = os.path.split(fp)
	return path, File
