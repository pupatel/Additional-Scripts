## -*- coding: utf-8 -*-
#"""


#Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
#Date created: 12/08/2015
##This script generates heatmapsclusters


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd




df=pd.read_csv('LATEST_V5.Parth.txt',sep='\t',index_col=0)
#df= df.transpose()



#cmap = sns.cubehelix_palette(as_cmap=True, dark=0.55,light=0.95) #this is the main pallete for figure 1
#cmap = sns.cubehelix_palette(as_cmap=True,dark=0.55,light=0.95)#, dark=0.55,light=0.95) #this is the main pallete for figure 1
#cmap=sns.color_palette("OrRd9")
#fig1=sns.clustermap(df,standard_scale=0,row_cluster=True,col_cluster=True,linewidths=.05,yticklabels=True,figsize=(15, 15),cmap=cmap) #this is the main pallete for figure 1

hmcol = ["#FFF9F0","#FEEDD3","#FDDDB1","#FDC99D","#FDA47A","#F2846D","#F2846D","#E94849"]#"#C23333"]#,"#993333"] # original
#hmcol = ["#ffffff","#fbe576","#c06e36","#9a2651"]
cmap = sns.blend_palette(hmcol,n_colors=1024,as_cmap=True)

#heahmapfunction
sns.set(font_scale=0.3) #0.3 for normal use
#fig1=sns.clustermap(df,cbar=True,standard_scale=0,row_cluster=True,col_cluster=False,linewidths=.05,yticklabels=True,figsize=(8.27, 11.69),cmap=cmap) # for landscape figsize=(11.69,8.27) or you can also multiply by it's ratio if you want bigger figure.
#fig1.cax.set_visible(False) #remove color bar
# for the paper add this fig1=sns.clustermap(df,standard_scale=0,row_cluster=True,col_cluster=False,linewidths=.05,square=True,yticklabels=True,figsize=(8.27, 11.69),cmap=cmap) 




###################Generated familyCounts_new_LATEST_V4_MORE_3_species.PDF on 10/10/2016
#
fig3=sns.clustermap(df,cbar=True,standard_scale=0,row_cluster=True,col_cluster=False,linewidths=.01,yticklabels=True,square=False,figsize=(15, 30),cmap=cmap) # for 24-mers # used on 4/5/2016
fig3.cax.set_visible(False) #remove color bar

hm = fig3.ax_heatmap.get_position()
plt.setp(fig3.ax_heatmap.yaxis.get_majorticklabels(), fontsize=1.5)
fig3.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.2, hm.height])    
col = fig3.ax_col_dendrogram.get_position()
fig3.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height*0.5])
#plt.show()


fig3.savefig('familyCounts.pdf',orientation='portrait',dpi=600)

#############################################################################################




print ("Done")











