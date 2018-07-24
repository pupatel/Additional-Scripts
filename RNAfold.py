
#!/usr/bin/python
#!/usr/bin/python
# coding=utf-8

#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 7/21/2015

##This script generates 2nd-ary structure from miRNA precursor using RNAFold and annotate  (highlight the miRNA) using RNAplot from Vienna RNA package.




import sys
import os
import subprocess





    
def split_at(s, c, n):
    words = s.split(c)
    return c.join(words[:n]), c.join(words[n:])
    
    
    
def RNAFOLD(NAME,LONG,SHORT):
#    pipe = subprocess.Popen(["perl","HPAnnotate.pl",str(LONG),str(SHORT),str(NAME)],stdout=subprocess.PIPE) # RUN PERL SCRIPT
    retcode2 = subprocess.call(["perl","HPAnnotate_test.pl",str(NAME),str(LONG),str(SHORT)]) # RUN PERL SCRIPT



if __name__ == '__main__':#calls the main function

    LONG_SEQ_FILE = str(sys.argv[1]) #"Aspa_mapped_genomic_seq.fa" 
    SHORT_SEQ_FILE = str(sys.argv[2]) #"SHORT.txt"
    
    fl_in = open(LONG_SEQ_FILE,'r')
    fs_in = open(SHORT_SEQ_FILE,'r') 
    lines = fs_in .readlines() #read SHORT SEQ. FILE INTO ARRAY
    data = fl_in.read().split('>')


    
  
    count=0 #loop counter
 
    for i in data[1:]:
        
        block = i.split('\n')
        name = block[0].split(' ')[0]        
        name = split_at(name,'_',2)
        #print name[0] 
        Long_SEQ = block[1]
        # print Long_SEQ
        Short_SEQ = lines[count].rstrip() #remove newline from python#fs_in.read().split('\n')[0]
        # print Short_SEQ
        # print name[0]
        RNAFOLD(name[0].rstrip(),Long_SEQ,SHORT_SEQ_FILE)
        count+=1
        
        print "Done..\n"
        
    
    
    
