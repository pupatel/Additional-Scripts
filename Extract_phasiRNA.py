#!/usr/bin/python
# coding=utf-8

#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 3/30/2015 | Date updated: 

##This script calculates number of phasiRNAs per phaseLoci with register length  of (by default 10 on both strands=20, so reg(20)x21=442)and extract them.

import MySQLdb
import numpy as np
#import scipy as sp
import matplotlib as mpl
import sys
## agg backend is used to create plot as a .png file
mpl.use('agg')
from pylab import *
import matplotlib.pyplot as plt
import random
from operator import itemgetter
from collections import defaultdict


    
# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
database_name="XXXXXXX"
user_name="XXXXX"
user_password="XXXXXX"
server="raichu.ddpsc.org" #"localhost"
nthPhas_length= 24 # or 21
ntheREGISTER=10  
nthshift=0 # deafualt 0 is recommended. mismatched reads (-1/+1), you should use (start+/- 1) to query against our dbs.
nthLIMIT= 50# num. of phasiRNAs to be extracted



def main(inputfile,out1):
      
    
    
    #connect to MySQL database
    db = MySQLdb.connect(host= server, user= user_name, passwd=user_password,db=database_name )
    cursor = db.cursor()
    
    inp_file_name = inputfile#str(sys.argv[1]) #user defined name for the input file
    out_file_name=out1
    out= open(out_file_name,'w');
    out.write('tag \t abundance \t strand \t cycle \t PHAS_ID \n')
    
    
#    phase_ID=""
    chr_ID=0
    start=0
    end=0
    
    
    
    INPUT= open(inp_file_name).readlines()
    #firstLine = INPUT.pop(0) #removes the first line
    num_lines = sum(1 for line in open(inp_file_name))
    print "num_of_PHASloci = "+ str(num_lines)
    
 
    phasiRNA=[]
    sorted_phasiRNA=[]
    
    line_index=0
    for lines in INPUT: #reading each line from the file - reading each PHAS locus 
        block = lines.split('\n')
        print block[0]
        phase_ID = block[0].split('\t')[0]
        chr_ID= int(block[0].split('\t')[1])
        start= int(block[0].split('\t')[2])+nthshift
        end= int(block[0].split('\t')[3])
#        PHASLoci_length=end-start
        
     
        index = ntheREGISTER
        register_increment=nthPhas_length
       
        watson_reg_beg=start
        watson_reg_end=0
        creek_reg_beg=0
        creek_reg_end=0
        
        for i in range (0,index):
            #print i
            #  extract phasiRNAs per register per PHASlocus on WATSON STRAND (UPPER)
            watson_reg_beg=start+(register_increment*i)
            watson_reg_end=watson_reg_beg+register_increment-1            
            # cursor.execute("select tag,norm_sum,strand FROM tag_position WHERE chr_id=%i and strand ='w' and (position between %i and %i) and length(tag)=%i order by norm_sum desc limit %i" %  (chr_ID,watson_reg_beg,watson_reg_end,nthPhas_length,nthLIMIT))
            cursor.execute("select tag,norm_sum,strand FROM tag_position WHERE chr_id=%i and strand ='w' and (position between %i and %i) and length(tag)=%i order by norm_sum desc" %  (chr_ID,watson_reg_beg,watson_reg_end,nthPhas_length))
	    rows=cursor.fetchall()
            #print rows
            for row in rows:
                tag,abundance,strand=row
                abundance=int(abundance)
                #print "w",tag,abundance               
                phasiRNA.append([tag,abundance,strand,i,phase_ID])
            
         
        
            
            
            
            #  extract phasiRNAs per register per PHASlocus on CREEK STRAND (LOWER)
            creek_reg_beg=(start+register_increment*i)+(register_increment-3) # -3 is for the 2-nt overhang
            creek_reg_end=creek_reg_beg-register_increment+1
            # cursor.execute("select tag,norm_sum,strand FROM tag_position WHERE chr_id=%i and strand ='c' and (position between %i and %i) and length(tag)=%i order by norm_sum desc limit %i" %  (chr_ID,creek_reg_end,creek_reg_beg,nthPhas_length,nthLIMIT))
            cursor.execute("select tag,norm_sum,strand FROM tag_position WHERE chr_id=%i and strand ='c' and (position between %i and %i) and length(tag)=%i order by norm_sum desc" %  (chr_ID,creek_reg_end,creek_reg_beg,nthPhas_length))
	    rows_1=cursor.fetchall()
            for row_1 in rows_1:
                tag_1,abundance_1,strand=row_1
                #print "c",tag_1,abundance_1
                abundance_1=int(abundance_1)
                phasiRNA.append([tag_1,abundance_1,strand,i,phase_ID])
            
     
            #for i in range(0,len(phasiRNA)):
                #print phasiRNA[i]
                #print "\n"
            
        line_index+=1
        
        
    temp={}
    for v in sorted(phasiRNA, key=lambda phasiRNA: phasiRNA[1]): # where L is your list
        temp[v[0]] = v
        print temp[v[0]]
    sorted_phasiRNA = temp.values()
    sorted_phasiRNA=sorted(sorted_phasiRNA,key=itemgetter(1),reverse=True)
    
    for i in range(0,nthLIMIT*num_lines):
            tag,abun,strand,register,phase_ID=sorted_phasiRNA[i]
            out.write('%s \t %d \t %s \t %d \t %s \n' % (tag,abun,strand,register,phase_ID))
        
    out.close()
 
   
        
if __name__ == '__main__':#calls the main function
    inp_file_name = "21_phaseLoci_Coordinates.txt" #user defined name for the input file with name,chromosome,start,end TAB seperated. (e.g. 21PHAS_NO8   1   159980586   159980986)
    out_file_name_1="24_Extracted_Rui.out" # Output file name   
    main(inp_file_name,out_file_name_1)
    print ("Done!\n")
    sys.exit()
