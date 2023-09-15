#!/usr/bin/env python
# Author: Ike Sanderson <ikes@uoregon.edu>

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
'''

__version__ = "0.3.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = None
RNA_bases = None

def reverse_complement(index2) -> str:
    '''reads sequence of nuceotides and returns their complementary sequence'''
    complementarydict={"A":"T","T":"A","C":"G","G":"C","N":"N"}
    revcomp=""
    index2=index2[::-1]
    for base in index2:
        revcomp+=(complementarydict[base])
    #print(f'{revcomp=}')
    return revcomp

def calc_median(sorted_list: list) -> float: 
    position1 = 0
    position2 = 0
    median = 0
    
    if len(sorted_list) % 2 == 1: #odd
        position1 = (len(sorted_list) // 2) #if odd, floor div length of list in half to ret mid point
        median = sorted_list[position1]  

    if len(sorted_list) % 2 == 0: #even    
        position1 = (len(sorted_list) // 2) #pos above med
        position2 = position1 - 1 #to get second, lower pos
        median = (sorted_list[position1] + sorted_list[position2]) / 2
    return median

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred 33 score"""
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """calculate the average quality score of the whole phred string"""
    qual_score=0
    i=0
    for value in (phred_score):
        score=convert_phred(value)
        qual_score=score+qual_score
        #print(f"{i}: {value} - {score} - {qual_score}")
        i+=1
    return qual_score/i

def validate_base_seq(sequence,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    #need to come back and fix this - go check the Python 2 day notes(?)
    seq = sequence.upper()
    DNAbases = set('ATGCNatcgn')
    RNAbases = set('AUGCNaucgn')
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")

def init_list(init_list: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list called init_list
    with 101 values of 0.0.'''

    for i in range(0,101):
        init_list.append(0.0)
    #else: 
        #print(len(init_list))#this is not needed for assignment but give me length info
    return init_list

def gc_content(DNA):
    '''how much G-C is there?'''
    DNA = DNA.upper()  #this is just in case you get lower case bases in the string
    CsnGs = DNA.count("C") + DNA.count("G") #count Gs and Cs
    return CsnGs / len (DNA) #calc percent

def oneline_fasta(filename):
    '''make a fasta file into a file with a header and sequence without endline chars'''
  # reads file NOTE: must write filename in argparse script NOTE and writes to outfile, pulling out the header while joining
    #all lines of the sequence
    header = ""
    sequence = ""
    with open(filename,'r') as ph, open("tempfile.fa",'w') as wh:
        wh.write(ph.readline())
        for line in ph:
            #line=line.strip('\n')
            if line.startswith('>'):
                if sequence == '':           
                    wh.write(f'\n{line}')
            else:    
                line = line.strip()
                wh.write(line) 

if __name__ == "__main__":
    #tests for reverse complement
    assert reverse_complement("A") == "T"
    assert reverse_complement("T") == "A"
    assert reverse_complement("C") == "G"
    #tests for validate_base_seq
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")
    # Leslie has already populated some tests for convert_phred - yay Leslie!
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")
    #tests for GC content
    assert gc_content("GCGCGC") == 1 #decimal/float of GC content
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    print("correctly calculated GC content") 