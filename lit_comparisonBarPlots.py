# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 21:20:09 2019

@author: twsle
"""

get_ipython().magic('reset -sf')


import matplotlib.pyplot as plt

import os
import networkx as nx
import pandas as pd
import numpy as np
import rdkit as rd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit import DataStructs
from matplotlib.ticker import FuncFormatter

from rdkit.Chem.Draw import DrawingOptions
from PIL import Image
from PIL import ImageDraw
from PIL import ImageDraw
from PIL import ImageFont
from PIL import ImageFile

os.chdir('C:\ResearchWorkingDirectory\CodeReferenceFiles_Python')

from CustomFunctions_PAHteam import FindSourceSink
from CustomFunctions_PAHteam import hierarchy_posMS


##Categories
#df = pd.read_excel(r"C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\SimiliarityCompounds.xlsx", sheet_name = 'HighThroughput')
#
#LogP = pd.read_excel(r"C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\AnthraceneRawQFnodes.xlsx", sheet_name = 'Sheet1')
#
#knownNodesDF = pd.read_excel(r"C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\SimiliarityCompounds.xlsx", sheet_name = 'Anthracene')
#
#listNodes = list(knownNodesDF['Perfect_Match'])
#knownNodes = [x for x in listNodes if str(x) != 'nan']
#
#LogP.drop(columns=['OutDeg'], inplace = True)
#LogP.drop(columns=['InDeg'], inplace = True)
#
#LogP2 = LogP.loc[LogP['nodeTransferWeights'] >= 0.01]
#LogP2.reset_index(inplace = True, drop = True)
##nodeList is the nodes predicted by the high Throughput method
#nodeList = list(LogP2['SMILES'])
#
##***********************
#finalNodes = []
#for n in nodeList:
#    m = Chem.MolFromSmiles(n)
#    #only include the edges if there's a ring
#    if (Chem.Lipinski.NumAromaticRings(m) != 0): 
#        finalNodes.append(n)
##***************************************
#
#df_raw = pd.read_excel(r"C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\AnthraceneRawWeight.xlsx", sheet_name = 'Sheet1')
#
##need to make sure we use the weights here
#G_raw = nx.from_pandas_edgelist(df_raw, source ='From', target = 'To', edge_attr=True, create_using=nx.DiGraph())
#
#newDF = df_raw.iloc[0:0]
#allNodes = []
#
##first nodelist, need this to generate Pathnodes
#for n in nodeList: 
#    m = Chem.MolFromSmiles(n)
#    #only include the edges if there's a ring
#    if (Chem.Lipinski.NumAromaticRings(m) != 0): 
#        #*****************************
#        outEdges = list(G_raw.out_edges(n))
#        tempc = []
#
#        for o in outEdges: 
#            if (o[0] in nodeList) &  (o[1] in nodeList):
#                tempc = df_raw.loc[(df_raw['From'] == o[0]) & (df_raw['To'] == o[1])]
#                newDF = newDF.append(tempc)
#                tempc = []
#            else:
#                print()
#        #*******************************
#
#
#G_new = nx.from_pandas_edgelist(newDF, source ='From', target = 'To', edge_attr=True, create_using=nx.DiGraph())
#
#[SourceN, SinkN] = FindSourceSink(G_new)
##drop the ghosts from Anthracene - these are the two compounds that are not part of the pathNodes 
#newDF.drop([45878,48358], inplace = True)
#
##dump the satillite nodes that are not part of path from the source:
#
##so now we can put together a list of the PathNodes
#pathNodes = list(set(newDF['From']) |set(newDF['To']) )
#
#
##*************************
#iterator = 0
#litNodeClassification = pd.DataFrame(columns=['SMILES', 'Classification'])
#needmatches = 0
#for k in knownNodes:
#    if k in nodeList:
#        litNodeClassification.loc[iterator] = [k,'perfectMatch']
#        iterator +=1
#    #check if the literature node matches the carbon backbone of anything in the model
#    elif k in list(LogP['SMILES']):
#        litNodeClassification.loc[iterator] = [k,'LowTPMatch']
#        iterator +=1
#    else:
#        litNodeClassification.loc[iterator] = [k,'NoMatch']
#        iterator +=1
#        #if there's no match, do through the PathNodes and dry to find the closest similiarty and then visually check the carbon backbone
#        m1 = Chem.MolFromSmiles(k)
#        fps1 = FingerprintMols.FingerprintMol(m1)
#        
#        #check the best possible match, and don't use path nodes, use the original list with the isolated ghosts in it
#        alldiceSims = pd.DataFrame(columns=['SMILES', 'Dice'])
#        iterator2 = 0
#        for n in nodeList:
#            #check the similitary
#            m2 = Chem.MolFromSmiles(n)
#            fps2 = FingerprintMols.FingerprintMol(m2)
#            diceSim = DataStructs.DiceSimilarity(fps1, fps2)
##            diceSim = DataStructs.FingerprintSimilarity(fps1, fps2)
#            alldiceSims.loc[iterator2] = [n,diceSim]
#            iterator2 += 1
#        alldiceSims.sort_values(by = ['Dice'], ascending = False, inplace = True)
#        alldiceSims.reset_index(drop = True, inplace = True)
#        print(k, ' best match is ', alldiceSims.loc[:3])
##        for i in range (0,4):
##            mOut = Chem.MolFromSmiles(alldiceSims.loc[i]['SMILES'])
##        mOut
#        needmatches +=1
##*******************************

#exact literature match with smiles
perfectMatches = []
#match against pathways in 
pathMatches = []
partialMatches = []
#missed a value from the literature, false negative
litNotModel = []
#predicted an extra value, false positive
modelNotLit = []


#'Acenapthene', 'Anthracene','Fluorene',  'Phenanthrene'
N = 4
plt.rcdefaults()
plt.rcParams["hatch.linewidth"] = 4
compoundTotals = np.array([16,47,46,84])

noMatchRaw = np.array([2,13,8,25])
noMatch =  np.divide(noMatchRaw,compoundTotals)

noMatchMultRaw =  np.array([0,6,4,13]) 
noMatchMult =  np.divide(noMatchMultRaw,compoundTotals)

partialLowRaw = np.array([0,3,5,4])
partialLow = np.divide(partialLowRaw,compoundTotals)

partialMatchRaw = np.array([4,11,11,25])
partialMatch = np.divide(partialMatchRaw,compoundTotals)

lowThroughputRaw = np.array([8,3,8,3])
lowThroughput = np.divide(lowThroughputRaw,compoundTotals)

perfectMatchRaw = np.array([2,11,10,14])
perfectMatch = np.divide(perfectMatchRaw,compoundTotals)

ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence


# red #CF3B2A
# blue #33A3D4
# green #499562
# grey #B1B1B2
# yellow #FCD238



#Set the grey hatch 


            
plt.rcParams['hatch.color'] = '#B1B1B2'      #  grey hatch       
#Perfect Match
p1 = plt.bar(ind, perfectMatch, label = 'Perfect Match',  color ='#499562') 
             
#Partial Match
p2 = plt.bar(ind, partialMatch, label = 'Partial Match',  bottom=perfectMatch, color ='#499562', hatch = r"//") #green with grey hatch 

    
                      
#Low Throughput  Perfect
p3 = plt.bar(ind, lowThroughput, label = 'Low Throughput Match', bottom=perfectMatch+partialMatch, color ='#33A3D4') #blue without grey hatch             
   
#Low Throughput Partial
p4 = plt.bar(ind, partialLow, label = 'Partial Low Throughput Match', bottom=perfectMatch+partialMatch+lowThroughput, color ='#33A3D4', hatch = r"//")  #blue with grey hatch
                         
#No Match Multiple Paper 
p5 = plt.bar(ind, noMatchMult, label = 'No Match Multiple Papers', bottom = perfectMatch+partialMatch+lowThroughput+partialLow,  color ='#CF3B2A')           #red color without hatch
             
plt.rcParams['hatch.color'] = 'black'                
#No Match Single Paper 
p6 = plt.bar(ind, noMatch, label = 'No Match Single Paper', bottom = perfectMatch+partialMatch+lowThroughput+partialLow + noMatchMult, facecolor='#CF3B2A', hatch = r"//")   # red color with hatch

#******************      
ax = plt.subplot(111)  
#ax.set_facecolor('white')
       
#ax = plt.gca()

box = ax.get_position()
ax.set_position([box.x0, 0.25, box.width * 0.65, box.height*.9])           
             
handles,labels = ax.get_legend_handles_labels()

handles = [handles[5], handles[4], handles[3], handles[2], handles[1], handles[0]]
labels = [labels[5], labels[4], labels[3], labels[2], labels[1],labels[0]]

#plt.legend(loc="best") 
plt.legend(handles,labels, loc='center left', bbox_to_anchor=(1, 0.5), facecolor = 'white',fontsize=14)  
#plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
#          ncol=3, fancybox=True, shadow=True)
#        
plt.xticks(ind, ('Acenapthene', 'Anthracene','Fluorene',  'Phenanthrene'), rotation = -45, fontsize=14)


#This deal is for the text on the plot
iterator = 0 
for r1, r2, r3, r4, r5, r6 in zip(p1, p2, p3, p4, p5, p6):
    h1 = r1.get_height()
    h2 = r2.get_height()
    h3 = r3.get_height()
    h4 = r4.get_height()
    h5 = r5.get_height()
    h6 = r6.get_height()
    
    
    #Perfect Match
    if perfectMatchRaw[iterator]>0:    
        plt.text(r1.get_x() + r1.get_width() / 2., h1/ 2., round(perfectMatchRaw[iterator],2), ha="center", va="center", color="white", fontsize=14, fontweight="bold")
    
    
    #Partial Match
    if partialMatchRaw[iterator]>0:    
        plt.text(r2.get_x() + r2.get_width() / 2., h1 + h2/ 2., round(partialMatchRaw[iterator],2), ha="center", va="center", color="white", fontsize=14, fontweight="bold")
    
    #LowThroughput Perfect Match
    if lowThroughputRaw[iterator]>0:     
        plt.text(r3.get_x() + r3.get_width() / 2., h1 + h2 + h3 / 2., round(lowThroughputRaw[iterator],2), ha="center", va="center", color="white", fontsize=14, fontweight="bold")
    
    #Low Throughput Partial
    if partialLowRaw[iterator]>0:         
        plt.text(r4.get_x() + r4.get_width() / 2., h1 + h2 + h3 +h4 / 2., round(partialLowRaw[iterator],2), ha="center", va="center", color="white", fontsize=14, fontweight="bold")    
    
    #No Match Multiple
    if round(noMatchMultRaw[iterator])>0:        
        plt.text(r5.get_x() + r5.get_width() / 2., h1 + h2 + h3 + h4 +h5 / 2., round(noMatchMultRaw[iterator],2), ha="center", va="center", color="white", fontsize=14, fontweight="bold")    
    
    #No Match Single
    if noMatchRaw[iterator]>0:
        plt.text( r6.get_x() + r6.get_width() / 2., h1 + h2 + h3 + h4 +h5 +h6/ 2., round(noMatchRaw[iterator],2), ha="center", va="center", color="white", fontsize=16, fontweight="bold")
    
    iterator += 1
    print(iterator)

ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
ax.tick_params(axis="y", labelsize=16)
ax.tick_params(axis="x", labelsize=16)

fig = plt.gcf()
cur_axes = plt.gca()

#fig.tight_layout()
fig.subplots_adjust(bottom = 0.3,top = 0.9, right = 0.6)


fig.set_size_inches(10, 6)

#fig.set_size_inches(2, 5)

for _, spine in ax.spines.items():
    spine.set_visible(True)
    spine.set_linewidth(1)
    spine.set_edgecolor('k')

plt.show()

os.chdir('C:\ResearchWorkingDirectory\TempImages')

fig.savefig('BarPlot2.png', format='png', dpi=2400)
#
import numpy as np
Image.MAX_IMAGE_PIXELS = None
image=Image.open('BarPlot2.png')
image.load()

#image_data = np.asarray(image)
#blankSpace = np.full((1000,20000,4), 0, dtype = 'uint8')
#
#image_data_final = np.vstack( (blankSpace,image_data_new))
#
font = ImageFont.truetype("arial.ttf", 550)

ImageDraw.Draw(image).text((100, 200),  'Metabolites Identified in Empirical Literature' , (0, 0, 0) ,font = font)
#
##draw.text((0, 0),"Sample Text",(255,255,255),font=font)
#image.save('BarPlot2_final.png')
#image.close()
#
