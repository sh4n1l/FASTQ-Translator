#!/usr/bin/env python3
# Name: Shreyas Anil (shanil) and Shivani Rao (surao)
#Shivani submitted the paper

import sys
class FastqReader :
    """
    Parse through a FASTQ file, identifying and yielding the ID line, sequence line, and quality score line.
    Based on Dr. Bernick's FastAreader class, uses the same methods with small tweaks.
    """
    def __init__ (self, fname=''):
        '''constructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)
 
    def readFasta (self):
        '''Parse through each line of the file, and determine if it is an ID, Sequence, or quality score line, or none of the above.'''
        ID = ''
        sequence = ''
        quality = ''
        
        with self.doOpen() as fileH:
            
            ID = ''
            sequence = ''
            quality = ''
 
            # skip to first fastq header
            line = fileH.readline()
            while not line.startswith('@') :
                line = fileH.readline()
            ID = line[1:].rstrip()
            plus = False #setting a variable to help decide between the sequence and quality lines
 
            for line in fileH:
                if plus == True:
                    quality = line.rstrip()
                    plus = False
                    continue
                if line.startswith ('@'): #ID line, should be the start of a new block
                    yield ID,sequence,quality
                    ID = line[1:].rstrip()    
                    sequence = ''
                    quality = ''
                elif line.startswith('+'): #next line will be a quality line 
                    plus = True 
                    continue
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
 
                 
        yield ID,sequence,quality
       
    
import re
import math
class FastqTranslator() :
    """
    Create the FastQtranslator class, which has methods to translate between different scoring systems for FastQ files.
    Take in a FastQ sequence and quality score as input, with either phred+33, phred+64, phred+64 with B as the lowest score,
    or phred+64 with the Solexa scoring matrix.
    Output a translated quality score with either phred+33 or phred+64.
    Design plan:
        Use 4 methods. 3 to convert from input to phred+33, and then a dictionary to convert from phred+33 to phred+64.
        This requires an extra step for going from phred+64 with B and phred+64 Solexa to the generic phred+64, but since
        it uses a dictionary, the extra step is quick. This way, the code is more efficient, with 2 less methods.
    Assumptions:
        User knows and input what type of scoring system the input uses. 
        Upper bound on error probability does not go higher than 40.
    """
    def __init__(self):
        '''
        Create a FastQtranslator object.
        Initialize a dictionary for converting from Sanger+33 to Illumina+64.
        '''
        self.sangIllDict = {}
        for i in range(33, 74):
            sang = chr(i)
            ill = chr(i + 31)
            self.sangIllDict[sang] = ill
            
    def p64ToP33(self, thisQual):
        '''Take in a phred+64 quality score and convert it to phred+33.'''
        tranQual = ''
        for score in (thisQual):
            num = ord(score) #find decimal representation of that ASCII character
            char = chr(num - 31) #change to 
            tranQual = tranQual + char
            
        return tranQual
    
    def p64bToP33(self, thisQual, thisSeq):
        '''Take in a phred+64 score with B as the lowest score, and convert it to phred+33.'''
        tranQual = ''
        revQual = reversed(thisQual) #quality score backwards
        pattern = re.compile('B+?') #regular expression, searching for string of one or more 'B's' non-greedily
        myQual = ''

        if pattern.match(str(revQual)) == 'None': #if the last letter of the quality score is not 'B'
            pass
        else:
            for score in revQual:
                if score == 'B':
                    myQual = myQual + "@" #make a string of @s, add one '@' sign for each 'B'
                else:
                    break 

        newQual = ''
        myLength = len(thisQual) - len(myQual) #accounts for any unknown bases at the end of the string
        for num in range(0, myLength): 
            if thisQual[num] == 'B' and thisSeq[num] == 'N': #if the quality score should really be a 0(error probability of 1)
                newQual = newQual + '@' #replaces the B with an @
            else:
                newQual = newQual + thisQual[num]
        finalQual = newQual + myQual #concatenate string with B's replaced with the string of @'s from before
                
        finalAns = self.p64ToP33(finalQual) #get the p33 output
        
        return finalAns            
    
    def p64SolToP33(self, thisQual):
        '''Take in a phred+64 score using the Solexa formula, and convert it to phred+33.'''
        tranSeq = ''
        for score in (thisQual):
            num = (ord(score)-64) #find decimal number associated with that ASCII character, this equals error probability
            num2 = int((math.log((10 ** (num / 10) + 1), 10)) * 10) #formula to convert between Solexa and phred scoring matrices
            if num2 <= 1:
                tranSeq = tranSeq + chr((num2) + 32)
            else:
                tranSeq = tranSeq + chr((num2) + 33) #return the character value associated with the new error probability
        return tranSeq
            
    def p33ToP64(self, thisQual):
        '''Take in a phred+33 score, and convert it to phred+64, using the dictionary from the __init__ method.'''
        tranSeq = ''
        for score in (thisQual):
            ill = self.sangIllDict.get(score, ' ') #used get() because direct dictionary lookup would not work on terminal for some reason
            tranSeq = tranSeq + ill
        return tranSeq

class CommandLine() :
    '''
    Uses Dr. Bernick's CommandLine class as template.
    Handle the command line, usage and help requests, using argparse.
    Implement a standard command line argument parser with various argument options,
    a standard usage and help.
    6 arguments: Phred+33 input, phred+64 input, phred+64 solexa input, phred+64 with B offset input,
        phred+33 output (default), and phred+64 output
    '''
    def __init__(self, inOpts=None) : 
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Enter the type of input file, and it can be translated into either Sanger Phred+33 or Illumina Phred+64', epilog = 'Thanks for using our program', add_help = True, #default is True 
                                              prefix_chars = '-', usage = '%(prog)s [options] -option1[default] <input >output' )
        #self.parser.add_argument('inFile', action = 'store', help='input file name') 
        #self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-P33in', '--PHRED33input', action = 'store', nargs='?', const=True, default=False, help='Sanger input')
        self.parser.add_argument('-P64in', '--PHRED64input', action = 'store', nargs='?', const=True, default=False, help='newer Illumina input')
        self.parser.add_argument('-P64Bin', '--PHRED64Boffset', action = 'store', nargs='?', const=True, default=False, help='PHRED64 with B offset in quality values') 
        self.parser.add_argument('-P64SOLin', '--PHRED64SOLEXA', action = 'store', nargs='?', const=True, default=False, help='PHRED64 with SOLEXA interpretation of Q score') #allows multiple list options
        self.parser.add_argument('-P33out', '--PHRED33output', action = 'store', nargs='?', const=True, default=True, help='Sanger output')
        self.parser.add_argument('-P64out', '--PHRED64output', action = 'store', nargs='?', const=True, default=False, help='newer Illumina output')
                                 
        if inOpts is None :
            self.args = self.parser.parse_args() 
        else :
            self.args = self.parser.parse_args(inOpts)
########################################################################
# Main
# Here is the main program
########################################################################
def main(inCL=None): 
    '''
    Translate a FASTQ file using either Sanger Phred+33 or Illumina Phred+64. 
    Take the command line input and feed it into the FastqTranslator class.
    Determine the output based on the arguments made on the command line.
    '''
    
    if inCL is None:
        myCommandLine = CommandLine() 
    else :
        myCommandLine = CommandLine(inCL)        
    
    print (myCommandLine.args)
    #takes file as stdin 
    myReader = FastqReader()
    newQual = FastqTranslator()
    for ID, seq, qual in myReader.readFasta():
        seq = seq.upper() #clean up the sequence
        seq = seq.replace('*', 'N')
        seq = seq.replace('.', 'N')
        seq = seq.replace('n', 'N')
        seq = seq.replace(' ', 'N')
        if (myCommandLine.args.PHRED64input == True) and (myCommandLine.args.PHRED64output == False): #Illumina input, Sanger output
            myScore = newQual.p64ToP33(qual)
            print ('@' + ID + '\n' + seq + '\n+\n' + myScore)
        elif (myCommandLine.args.PHRED33input == True) and (myCommandLine.args.PHRED64output == True): #Sanger input, Illumina output
            myScore = newQual.p33ToP64(qual)
            print ('@' + ID + '\n' + seq + '\n+\n' + myScore)
        elif (myCommandLine.args.PHRED64Boffset == True) and (myCommandLine.args.PHRED64output == False): #special Illumina input, Sanger output
            myScore = newQual.p64bToP33(qual, seq)
            print ('@' + ID + '\n' + seq + '\n+\n' + myScore)  
        elif (myCommandLine.args.PHRED64Boffset == True) and (myCommandLine.args.PHRED64output == True): #special Illumina input, normal Illumina output
            firstScore = newQual.p64bToP33(qual, seq)
            myScore = newQual.p33ToP64(firstScore)
            print ('@' + ID + '\n' + seq + '\n+\n' + myScore) 
        elif (myCommandLine.args.PHRED64SOLEXA == True) and (myCommandLine.args.PHRED64output == False): #Solexa input, Sanger output
            myScore = newQual.p64SolToP33(qual)
            print ('@' + ID + '\n' + seq + '\n+\n' + myScore)
        elif (myCommandLine.args.PHRED64SOLEXA == True) and (myCommandLine.args.PHRED64output == True): #Solexa input, Illumina output
            firstScore = newQual.p64SolToP33(qual)
            myScore = newQual.p33ToP64(firstScore)
            print ('@' + ID + '\n' + seq + '\n+\n' + myScore)

if __name__ == "__main__":
    main()


