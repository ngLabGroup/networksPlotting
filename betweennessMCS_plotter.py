# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 12:57:57 2019

@author: twsle
"""
get_ipython().magic('reset -sf')

import os
import numpy as np

import pandas as pd

#PIL imports
from PIL import Image    
from PIL import ImageDraw
from PIL import ImageFont
from PIL import ImageFile


import matplotlib.pyplot as plt
import networkx as nx

from IPython.display import display # to display images

from rdkit import Chem
from rdkit.Chem import AllChem
from matplotlib import colors

from rdkit.Chem import Descriptors
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdFMCS
#from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import DrawingOptions
from rdkit.Chem.Draw import FingerprintEnv

import xlsxwriter

#set the directoring
os.chdir('C:\ResearchWorkingDirectory\CodeReferenceFiles_Python')

#Define clustering setup
#*******************************************************
def ClusterFps(fps,cutoff=0.2):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs
#*******************************************************
    
#required input file: QF file from the nodeThroughput script
    

#This file will not work with the example data as it requires at least several hundred nodes in order to develop MCS structures. However, the code is provided as reference. 
dfQF = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\networksPlotting\PhenQFNodes20191226-175705.xlsx", sheet_name = 'Sheet1')

#sort and take the top 1% of highest betweeness compounds
dfQF.sort_values(by ='Betweenness Centrality', ascending = False, inplace = True)
dfQF.reset_index(inplace = True, drop = True)
top1Percent = round(len(dfQF)/100)

topSmiles = dfQF[['SMILES', 'Betweenness Centrality']].head(top1Percent)



#build the molecule list
smilesList = list(topSmiles['SMILES'])
ms = [Chem.MolFromSmiles(s) for s in smilesList]

fps = [AllChem.GetMorganFingerprint(x,3,useFeatures=True) for x in ms]
clusters=ClusterFps(fps,cutoff=0.4)

#sort the clusters so the largest ones are on top
cluster2 = sorted(clusters, key=len, reverse = True) 
       
#Add the cluster identification to the rows
for i, row in topSmiles.iterrows():
    i3 = 0
    while i not in cluster2[i3]:
        i3 +=1
        pass
    
    topSmiles.at[i,'Cluster Assignment'] = i3

#Draing options for the display of the MCS structures at the top of the file
DrawingOptions.atomLabelFontSize = 120
DrawingOptions.dotsPerAngstrom = 120
DrawingOptions.bondLineWidth = 20.0
DrawingOptions.colorBonds   =  False
DrawingOptions.noCarbonSymbols  = True
DrawingOptions.defaultColor = (1,0.5,0)
DrawingOptions.dblBondOffset = 0.35

allMCS = []  
allMCSimg = []  
#this look give you the MCS of each cluster

#8 controls oxygen
DrawingOptions.elemDict= {0: (0.5, 0.5, 0.5), 1: (0.55, 0.55, 0.55), 7: (0, 0, 1), 8: (1,0.5,0), 9: (0.2, 0.8, 0.8), 15: (1, 0.5, 0), 16: (0.8, 0.8, 0), 17: (0, 0.8, 0), 35: (0.5, 0.3, 0.1), 53: (0.63, 0.12, 0.94)}

for i2 in range (0,10):
    mols = [ms[i] for i in cluster2[i2]]
    
    res = Chem.rdFMCS.FindMCS(mols,completeRingsOnly=True)
    
    ss = Chem.MolFromSmarts(res.smartsString)   
    testS = res.smartsString
    
    mcsp = Chem.MolFromSmarts(res.smartsString)
    Chem.rdmolops.Kekulize(mcsp)

    newS = Chem.MolToSmiles(mcsp,kekuleSmiles=True) 
    print(newS)
    
    #artificially kekulize the smile string so it's consistent with everything else. 
    a = newS.replace(":C:C:C:C:C:", "-C=C-C=C-C=")

    m = Chem.MolFromSmiles(a)
    
    allMCS.append(m)
    
    img =  Chem.Draw.MolToImage(m, size=(1200, 1200),fitImage=False, Kekulize=True)
    img = img.crop((0,0,1200,1000))
    allMCSimg.append(img)
#    display(img)
    
#    print (len(cluster2[i2]))
    
#img = Draw.MolsToGridImage(allMCS, molsPerRow=5, useSVG=False)
#display(img)



#**********************************************************
tempIMGS = []
compoundsinTop5 = []

DrawingOptions.atomLabelFontSize = 120
DrawingOptions.dotsPerAngstrom = 120
DrawingOptions.bondLineWidth = 10.0
DrawingOptions.defaultColor = (0,0,0)

#Set up the elemDict so that we get the colors right
DrawingOptions.elemDict= {0: (0.5, 0.5, 0.5), 1: (0.55, 0.55, 0.55), 7: (0, 0, 1), 8: (0, 0, 0), 9: (0.2, 0.8, 0.8), 15: (1, 0.5, 0), 16: (0.8, 0.8, 0), 17: (0, 0.8, 0), 35: (0.5, 0.3, 0.1), 53: (0.63, 0.12, 0.94)}

for whichCluster in range(0,5):
    pass
    print('i am in ', whichCluster)
#    whichCluster = 4
    mols = [ms[i] for i in cluster2[whichCluster]]
    
    res = Chem.rdFMCS.FindMCS(mols,completeRingsOnly=True)
    
    ss = Chem.MolFromSmarts(res.smartsString)   
    testS = res.smartsString
    
    mcsp = Chem.MolFromSmarts(res.smartsString)
    
#    img = Draw.MolsToGridImage([mcsp])

    #Second loop
#    hcolor = colors.to_rgb('orange')
    hcolor = (1,0.5,0)
    mols = []
    atomsets = []
    hcolors = []
    allImages2 = []
    
    compoundsinTop5.append(len(cluster2[whichCluster]))
    
    for c in cluster2[whichCluster]:
        mol = ms[c]
        mols.append(mol)
        matches = mol.GetSubstructMatches(mcsp)
        atomset = list(matches[0])
        atomsets.append(atomset)
        hcolors.append(hcolor)
        
        img = Chem.Draw.MolToImage(mol, size=(1200, 1200),fitImage=False, highlightAtoms=atomset, highlightColor=hcolor)
#        display(img)
        allImages2.append(img)
    
    #*******************
    
    left = 0
    top = 100
    right = 1200
    bottom = 1100
    
    allImages = [a.crop((left,top,right,bottom)) for a in allImages2[0:5]]

    imgs_comb2 = np.vstack( (np.asarray( i ) for i in allImages ) )
    
    if len(allImages2) < 5:
        #Fill with a blank image
        blankImage = np.full((1000,1200,4), 255, dtype = 'uint8')
        #stack another one on if necessary
        imgs_comb2 = np.vstack( (imgs_comb2, blankImage))

#destroy this later
#    for i in allImages:
#        pass
#      
##    print(len(allImages))
##    #if one of the images is the wrong size halve to add a dummy image
##    if len(allImages) > 5:
##        pass
       
    tempIMGS.append(imgs_comb2)
    #add the MCS at the top

# save that beautiful picture
ImageFile.MAXBLOCK = 2**20

#need to stack a side border also
blackSide = np.full((5000,15,4), 255, dtype = 'uint8')


imgs_final = np.hstack( (tempIMGS[0],tempIMGS[1],tempIMGS[2],tempIMGS[3],tempIMGS[4],blackSide) )


tempMCS = allMCSimg[0:5]

blackSide = np.full((1000,15,4), 255, dtype = 'uint8')

imgs_top = np.hstack((tempMCS[0],tempMCS[1],tempMCS[2],tempMCS[3],tempMCS[4],blackSide ) )

#add a little bit of blank space on the top for the text
blankBorder = np.hstack( (np.full((550,6000,4), 255, dtype = 'uint8'), np.full((550,15,4), 255, dtype = 'uint8') ) )

blackBorder = np.full((15,6015,4), 0, dtype = 'uint8')

imgs_final2 = np.vstack( (blankBorder,imgs_top,imgs_final))

#need a white buffer border also
whiteBorder = np.full((6550,50,4), 255, dtype = 'uint8')

imgs_final3 = np.hstack( (whiteBorder,imgs_final2 ))

imgs_comb = Image.fromarray( imgs_final3)
imgs_comb = imgs_comb.convert('RGB')

font = ImageFont.truetype("arial.ttf", 220)

# Option to add text as a title 
#ImageDraw.Draw(imgs_comb).text((50, 0),'D. Phenanthrene', (0, 0, 0) ,font = font)

#Test labels
ImageDraw.Draw(imgs_comb).text((50, 225),'Number of Compounds in Top 1% Betweenness: ' + str(top1Percent), (0, 0, 0) ,font = font)

ImageDraw.Draw(imgs_comb).text((50, 450),'Number of Compounds in Top 5 Structure Clusters: ' + str(np.sum(compoundsinTop5)), (0, 0, 0) ,font = font)

#write the number of compounds int he cluster
ImageDraw.Draw(imgs_comb).text((100, 700),  str(compoundsinTop5[0])  , (0, 0, 0) ,font = font)
ImageDraw.Draw(imgs_comb).text((1300, 700), str(compoundsinTop5[1]) , (0, 0, 0) ,font = font)
ImageDraw.Draw(imgs_comb).text((2400, 700), str(compoundsinTop5[2]) , (0, 0, 0) ,font = font)
ImageDraw.Draw(imgs_comb).text((3600, 700), str(compoundsinTop5[3]) , (0, 0, 0) ,font = font)
ImageDraw.Draw(imgs_comb).text((4800, 700), str(compoundsinTop5[4]) , (0, 0, 0) ,font = font)

#draw lines
ImageDraw.Draw(imgs_comb).line( [(0,1550),(6015,1550)], fill = (0,0,0), width = 15)
ImageDraw.Draw(imgs_comb).line( [(0,700),(6015,700)], fill = (0,0,0), width = 15)

##save off that beautiful image
imgs_comb.save( 'Phen BetweennessTest.jpg' , "JPEG", quality=100, optimize=False, progressive=False)    
