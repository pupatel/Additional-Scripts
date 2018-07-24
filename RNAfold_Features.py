
#!/usr/bin/python
#!/usr/bin/python
# coding=utf-8

#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 5/10/2016 Date modified:7/24/2018

##This script generates 2nd-ary structure from miRNA precursor using RNAFold and parse the 2nd strcutre and generate structural features (dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE,Freq,Div) using RNAplot from Vienna RNA package.




import sys
import os
import subprocess
import itertools,re
import matplotlib.pyplot as plt



#global variable
Struct_Feat={} # dictionary to hold all structural features
total_bp=[]
num_stems=[]
num_loops=[]

    
def split_at(s, c, n):
    words = s.split(c)
    return c.join(words[:n]), c.join(words[n:])
    
    
    
def RNAFOLD(Seq):

    os.environ['Seq']=Seq
    os.system("echo $Seq| RNAfold -d2 -p > RNAfold.out")

def PARSER(File,count):
    fl_in = open(File,'r')
    lines = fl_in .readlines() #read SHORT SEQ. FILE INTO ARRAY

    #Compute total number of base pairs link - Dot bracket notation (http://rna.tbi.univie.ac.at/help.html)
    dot_notation=lines[1].split(" ", 1)[0]
    total_bases= dot_notation.count('(')
    total_bp.append(total_bases) # To compute average total basepairs, we need this array. Added on 12/1/2016
    # print "Total_bases=",total_bases

    #Compute total number of stems: n_stems is the number of stems in the secondary structure. A stem is a structural motif containing more than three contiguous base pairs
    dot_notation_parse=[(k, sum(1 for _ in vs)) for k, vs in itertools.groupby(dot_notation)] # e.g. [('(', 9), ('.', 3), ('(', 6), ('.', 9), (')', 6), ('.', 8), ('(', 6), ('.', 7), (')', 6), ('.', 2), (')', 9)]
    n_stems=0
    for i in range (0,len(dot_notation_parse)):
        if dot_notation_parse[i][0]=='(': # find all list with "(" which means there is a paired base
            if dot_notation_parse[i][1]>=3: #look for more than three contiguous base pairs
                n_stems+=1
    num_stems.append(n_stems) # To compute average stems, we need this array. Added on 12/1/2016




    #Compute total number of loops: https://www.biostars.org/p/4300/
    loops = []
    n_loops=0 # number of loops
    p = re.compile('[(][.]{1,100}[)]')
    loopsIter = p.finditer(dot_notation)
    # print loopsIter
    for m in loopsIter:    
        loops.append(m.span())
    n_loops=len(loops)
    num_loops.append(n_loops)
    # print "loops",len(loops)



    # Compute dP = 'Normalized base-pairing propensity'
    dP=float(total_bases)/len(lines[0].strip())
    # print "dP = ", dP
    #parse MFE
    lines[1]=lines[1].split(" ", 1)
    mfe= lines[1][1].split("(",1)[1].split(")",1)[0]
    # print "mfe = ",mfe

    #compute dG=norm_mfe
    dG=float(mfe)/len(lines[0].strip())
    # print len(lines[0].strip())
    # print "dG = ", dG

    #compute MFEI1 - 'Index 1 based on the minimum free energy'
    seq=lines[0].strip()
    g = seq.count("G")
    c = seq.count("C")
    t = seq.count("U")
    a = seq.count("A")
    gcCount = g+c
    totalBaseCount = g+c+t+a    
    gcFraction = float(gcCount) / totalBaseCount

    ##Compute MEFE1 - 'Index 1 based on the minimum free energy'
    if(gcFraction>0):
        MFEI1= dG/(gcFraction*100)
        # print "MFEI1 = ",MFEI1
    else:
        MFEI1 = dG

    ##Compute MEFE2 - 'Index 2 based on the minimum free energy'
    if(n_stems>0):
        MFEI2= float(dG)/n_stems
        # print "MFEI2 = ",MFEI2,n_stems
    else:
        MFEI2 = 0 

        ##Compute MEFE2 - 'Index 2 based on the minimum free energy'
    if(n_loops>0):
        MFEI3= float(dG)/n_loops
        # print "MFEI3 = ",MFEI3,n_loops
    else:
        MFEI3 = 0 
    
    #Compute MEFE4 - 'Index 4 based on the minimum free energy'
    if(total_bases>0):
        MFEI4= float(mfe)/total_bases
        # print "MFEI4 = ",MFEI4
    else:
        MFEI4 = 0 

    #parse EFE
    lines[2]=lines[2].split(" ", 1)
    efe= lines[2][1].split("[",1)[1].split("]",1)[0]
    # print "efe = ",efe

    #compute NEFE
    NEFE=float(efe)/len(lines[0].strip())
    # print "NEFE = ",NEFE

    #parse FREQ
    Freq=lines[4].split("ensemble ")[1].split(";",1)[0]
    # print "freq = ",Freq

    #parse Diversity
    Div=lines[4].split("diversity ",1)[1].strip()
    # print "Div = ",Div
    
    #compute Difference
    Diff=(float(mfe)-float(efe))/len(lines[0].strip())
    # print Diff


    Struct_Feat[count]=[dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE,Freq,Div,Diff]    
    fl_in.close()
    # return dP,dG,MFEI1,MFEI4,NEFE,Freq,Div


if __name__ == '__main__':#calls the main function


    SHORT_SEQ_FILE= "Input_Sequence_File.txt" #files containing precursor sequences (a sequence per line).

    

    fs_in = open(SHORT_SEQ_FILE,'r')

    out_file_name="Strucutal_feat.txt"#Strucutal_feat_norway_spruce_21.txt" #"test_result.txt"
    out= open(out_file_name,'wb');
    out.write('dp \t dG \t MFEI1 \t MFEI2 \t MFEI3 \t MFEI4 \t NEFE \t Freq \t Div \t Diff\n')

    lines = fs_in .readlines() #read SHORT SEQ. FILE INTO ARRAY
       
    
  
    count=0 #loop counter
    for i in range(0,len(lines)):
        seq=lines[i]      
        RNAFOLD(seq)
        PARSER("RNAfold.out",count)   
        count+=1 
         
    for i in range(0,len(lines)):
        dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE,Freq,Div,Diff=Struct_Feat[i]
        # print Struct_Feat[i]
        out.write('%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n' % (dP,dG,MFEI1,MFEI2,MFEI3,MFEI4,NEFE,Freq,Div,Diff))
        
    out.close()

    print "Max number of stems = ",max(num_stems)
    print "Max number of total_bp = ", max(total_bp)
    print "Max number of loops = ",max(num_loops)

    print "Min number of stems = ",min(num_stems)
    print "Min number of total_bp = ", min(total_bp)
    print "Min number of loops = ",min(num_loops)

    print "Average number of stems = ",float(sum(num_stems)/len(num_stems))
    print "Average number of total_bp = ", float(sum(total_bp)/len(total_bp))
    print "Average number of loops = ",float(sum(num_loops)/len(num_loops))

    # plt.boxplot(num_stems)
    # plt.savefig('num_stems.png')
    # plt.boxplot(total_bp)
    # plt.savefig('total_bp.png')
    # plt.boxplot(num_loops)
    # plt.savefig('num_loops.png')

    # data=[num_stems, total_bp, num_loops]
    # labels = list('num_stemstotal_bpnum_loops')
    # plt.boxplot(data,labels=labels, showmeans=True, meanline=True)

    # plt.savefig('strcture.png')


 
    
print "Done..Strucutal_features Generates!"   
