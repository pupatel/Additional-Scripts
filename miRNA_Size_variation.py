# -*- coding: utf-8 -*-
"""
Created on Mon Jun 06 14:10:18 2016

@author: Admin
"""

import csv,os,re,sys
from collections import Counter
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict

#with open('Final_Report_ALL_Libs_MERGED_SM_All_plant_miRs_Duplicates_removed.csv') as csvfile:
#    readCSV = csv.reader(csvfile, delimiter=',')
#    for row in readCSV:
#        print row[0]
    
f= open('Filtered_21316_reads_=100_miR165_to_miR166_USE_this.csv','r') #open('filtered_3729.csv','r') #used 'Final_Report_ALL_Libs_MERGED_SM_All_plant_miRs_Duplicates_removed_LATEST.csv' for all analysis before 12/16/16 after Sandra askem to add new pienapple libraries.
o= open("Size_Distribution_intermediate_conserved_miRs_abundance_List.csv",'wb')#open("miRNAs_ALL_SPECIES_reza_filtered_3729_TEMP.csv",'wb')

#Conserved_list=['miR165','miR401','miR1432','miR158','miR5077','miR5083','miR5059','miR3630','miR6478','miR8155','miR8175','miR858','miR6173','miR477','miR403','miR5368','miR4376','miR894','miR5072','miR1511','miR3522','miR1510','miR4414','miR5225','miR157','miR7767','miR6300','miR5139','miR170','miR5078']

#highly conserved
#Highly_Conserved_list=['miR396','miR156','miR171','miR166','miR319','miR167','miR172','miR159','miR169','miR164','miR168','miR393','miR399','miR408','miR160','miR390','miR397','miR394','miR529','miR827','miR398']#['miR165','miR401','miR1432','miR158','miR5077','miR5083','miR5059','miR3630','miR6478','miR8155','miR8175','miR858','miR6173','miR477','miR403','miR5368','miR4376','miR894','miR5072','miR1511','miR3522','miR1510','miR4414','miR5225','miR157','miR7767','miR6300','miR5139','miR170','miR5078']

#intermediate conserved
Highly_Conserved_list=['miR162','miR528','miR395','miR444','miR535','miR2275','miR1507','miR2118','miR5179','miR530','miR894','miR1432','miR8175','miR6478','miR482','miR8155']



csv_f=csv.reader(f,delimiter=',')
csv_o=csv.writer(o,delimiter=',')

headers = next(csv_f, None) 
#print headers

miR_conservation__20mers_List={}
miR_conservation__21mers_List={}
miR_conservation__22mers_List={}
miRNA_20_mers=0
miRNA_21_mers=0
miRNA_22_mers=0



write_row="miRNA_name","20 nt","21 nt","22 nt"
csv_o.writerow(write_row)

for row in csv_f:
    
####### parse miRNA name##############
    
    miRNA_key=row[1].strip()
#    print miRNA_key
              
                        
                        
####### get miRNA family and sizes in that family ##############
    size=int(row[9].strip()) 
    abundance= int(row[27])
        
#    size
    if miRNA_key in Highly_Conserved_list:
        
        if (size==20):        
            if miRNA_key in miR_conservation__20mers_List: 
                new_abundance=int(miR_conservation__20mers_List[miRNA_key])+abundance
                miR_conservation__20mers_List[miRNA_key]=new_abundance
            else:
                miR_conservation__20mers_List[miRNA_key]=abundance
    #            miR_conservation__20mers_List[miRNA_key].append(Abundance)
        elif (size==21):        
            if miRNA_key in miR_conservation__21mers_List: 
                new_abundance=int(miR_conservation__21mers_List[miRNA_key])+abundance
                miR_conservation__21mers_List[miRNA_key]=new_abundance
            else:
                miR_conservation__21mers_List[miRNA_key]=abundance
    #            miR_conservation__20mers_List[miRNA_key].append(Abundance)
        else:        
            if miRNA_key in miR_conservation__22mers_List: 
                new_abundance=int(miR_conservation__22mers_List[miRNA_key])+abundance
                miR_conservation__22mers_List[miRNA_key]=new_abundance
            else:
                miR_conservation__22mers_List[miRNA_key]=abundance
    #            miR_conservation__20mers_List[miRNA_key].append(Abundance)
            
#print miR_conservation__20mers_List,miR_conservation__21mers_List ,miR_conservation__22mers_List 
#sys.exit()

#key_list=['Sorghum','Wheat','Sugarcane','Brachypodium','Rice','Riceglab','Setaria','Maize','Streptochaeta','Anomochloa','Pharus','Raddia','BambooEdulis','BambooMoso','Bromeliad','Pineapple','OilPalm','Cyperus','Coconut','DatePalm','Banana','Canna','Tradescantia','Prosthechea','Phalaenopsis','Daylily','Kniphofia','Liriope','Asparagus','Freycinetia','Lilium','Spirodela','Lemna','Colocasia','Echinodorus','Sagittaria','Zostera','Acorus','Arabidopsis','Nymphaea','Amborella']
#new_miR_conservation_List = OrderedDict((k, miR_conservation_List.get(k)) for k in key_list)


for miRNA_key in Highly_Conserved_list:
    
    if miRNA_key in miR_conservation__20mers_List:
        abundance_20_mers=int(miR_conservation__20mers_List[miRNA_key])
    else: 
        abundance_20_mers=0
        
    if miRNA_key in miR_conservation__21mers_List:
        abundance_21_mers=int(miR_conservation__21mers_List[miRNA_key])
    else:
        abundance_21_mers=0
        
    if miRNA_key in miR_conservation__22mers_List:
        abundance_22_mers=int(miR_conservation__22mers_List[miRNA_key])
    else:  
        abundance_22_mers=0
        
    total_miRNA_abundance=abundance_20_mers+abundance_21_mers+abundance_22_mers
    
    
#    print length_list,miR_conservation_List[key],key,miR_conservation_List[key].count('20'),miR_conservation_List[key].count('22'),miR_conservation_List[key].count('21')
    write_row=miRNA_key,float(abundance_20_mers)/total_miRNA_abundance*100,float(abundance_21_mers)/total_miRNA_abundance*100,float(abundance_22_mers)/total_miRNA_abundance*100
#    print write_row,total_miRNA_abundance
    csv_o.writerow(write_row)
#    sys.exit()

o.close()
f.close()

sns.set(rc={'figure.figsize':(8,7)})
sns.set_style("white")
df=pd.read_csv('Size_Distribution_intermediate_conserved_miRs_abundance_List.csv',sep=',')
df=df.sort_values(by=['22 nt'],ascending =[True]) #df=df.sort_values(by=['22 nt','21 nt'],ascending =[True,True])
Palette=['#FF5E5E','#405888','#45C57F']#'']
#x=df.groupby(level=[0]).sum().sort_values(ascending=False)
#print x
ax=df.plot(x=df['miRNA_name'],kind='bar', stacked=True, alpha=0.9, color=Palette,legend=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
#df.reindex(index=x.index, level=0).unstack().plot.bar(stacked=True, cmap=Palette)
#For 3' 3_prime_all_species_conservation_List only since C is most baised at 3' end in all species, we want to show c first 
#my_df = df[["Organism_name","C","T","A","G"]] 
#Palette=['#FED700','#004D64','#19C878','#FE6347']#['#19C878','#FED700','#FE6347','#004D64']#'']
#    
#ax=my_df.plot(x=my_df['Organism_name'],kind='bar', stacked=True, alpha=0.6, color=Palette,legend=False)

    
ax.figure.savefig('Size_Distribution_intermediate_conserved_miRs_List_7_9_18.png',dpi=600)


print "Done"

    
   
 