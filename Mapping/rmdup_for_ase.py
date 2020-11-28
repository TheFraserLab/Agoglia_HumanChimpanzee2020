#!/usr/bin/env python3
'''
FUNCTION: Remove duplicates from alignment file

NOTES:
- Files must end with .bam or .sam and must be in the specified format.  (The program interprets all .bam endings to be BAM files, etc)
- Header files must be attached to the file. 
- Input file should be sorted by samtools.
- EDIT: To run on sherlock server, subprocess.Popen() must have a single string as the command argument and must set shell=True.
    (R. Agoglia, 07162016)

TO DO:
- Allow for piping in data. 
- Add single-end rmdup.
- Add testing for paired-end rmdup using unittest
- Add check for sorted
- Add check for samtools0.1.19
- Add ability to put back in the original quality score.  
'''

from __future__ import print_function
import random, sys, os, subprocess, argparse

__author__="Ryo Kita"
__data__="11.8.2015"
def main():

    args = run_argparse()
    rmdupFromAlignment(args.input, args.output, args.paired, args.formatIn, args.formatOut, args.samtoolsPath)

#########################
def run_argparse():
    parser = argparse.ArgumentParser(description="Remove Duplicates from Paired-End Data")
    parser.add_argument('-i', '--input', required=True,  help="Input file: It must be in BAM or SAM format and end with .bam or .sam")
    parser.add_argument('-o', '--output', required=True,  help="Output file: It must end in .bam or .sam.")
    parser.add_argument('-p', '--paired', dest='paired', action='store_true', help="Include -p for paired end")
    parser.add_argument('-fi', '--formatIn',  help="Format of the file.  Must be either 'bam' or 'sam'. Automatically detected if missing.")
    parser.add_argument('-fo', '--formatOut',  help="Format of the file.  Must be either 'bam' or 'sam'. Automatically detected if missing.")
    parser.add_argument('-s', '--samtoolsPath', default= "/home/rkita/home/extProgs/samtools-0.1.19/samtools"   , help="Path to samtools0.19.  Note, the newest samtools will not work.")
    args = parser.parse_args()

    if args.formatIn==None:
        ending = args.input.split(".")[-1]
        if ending not in ["bam", "sam"]:
            print("ERROR: No format detected (-fi). Specify -fi or end the input filename with '.bam' or'.sam'")
            return 0
        else: args.formatIn=ending
    
    if args.formatOut==None:
        if args.output==None or  args.output.split(".")[-1] not in ["bam", "sam"]:
            print("ERROR: No format detected (-fo). Specify -fo or end the output filename with '.bam' or '.sam'")
            return 0
        else: args.formatOut=args.output.split(".")[-1] 
    return args

#########################
def rmdupFromAlignment(inputFile, outputFile, paired, formatIn, formatOut, samtoolsPath):

    print("Running Remove Duplicates From Ase")
    print("OPTIONS BEING USED: ")
    print("Input File ", inputFile)
    print("Output File ", outputFile)
    print("Paired End? ", paired)
    print("Input File Type: ", formatIn)
    print("Output File Type: ", formatOut)
    print("Samtools0.19 Path: ", samtoolsPath)
    print("\n\n")

    #STEP 1: Randomize the quality
    tmpFileName = inputFile + ".tmp.bam"
    if formatIn=="bam":
        viewSam = subprocess.Popen(samtoolsPath + " view -h " + inputFile, stdout=subprocess.PIPE, bufsize=1, shell=True)
    else: 
        viewSam = subprocess.Popen("cat " + inputFile, stdout=subprocess.PIPE, bufsize=1, shell=True)
    convertBam = subprocess.Popen(samtoolsPath + " view -uhS -o " + tmpFileName + " -", stdin=subprocess.PIPE, bufsize=1, shell=True)
    for line in iter(viewSam.stdout.readline, b''):
        convertBam.stdin.write(randomize_quality(line.decode()).encode())
    viewSam.stdout.close()
    viewSam.wait()
    convertBam.stdin.close()
    convertBam.wait()

    #STEP 2: Remove duplicates. 
    if formatOut=="bam":
        rmdup = subprocess.call(samtoolsPath + " rmdup " + tmpFileName + " " + outputFile, shell=True)
    else: #If SAM output requested, then convert to sam file at the same time
        fileObj = open(outputFile, "w")
        rmdup = subprocess.Popen(samtoolsPath + " rmdup " + tmpFileName + " -" , stdout = subprocess.PIPE, shell=True)
        convertToSam =  subprocess.Popen(samtoolsPath + " view -h -", stdin = rmdup.stdout, stdout = fileObj, shell=True)
        rmdup.wait()
        convertToSam.wait()

    #STEP 3: Clean up
    os.system("rm " +  tmpFileName)


#########################
def randomize_quality(line, seed=False):
    '''Randomize the quality of a sam formatted line'''
    if seed:         
        #set seed to function input
        random.seed(seed)
    if line[0]!="@": #ignore, but return if header
        #convert the 4th index quality score.
        lineSplit = line.split("\t")
        lineSplit[4] = str(random.randint(235,255))
        line = "\t".join(lineSplit)
    return line

#########################
if __name__=="__main__":
    main()
#########################

