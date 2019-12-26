# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 21:20:09 2019

@author: twsle
"""

get_ipython().magic('reset -sf')



import os
import networkx as nx
import pandas as pd
import numpy as np
import rdkit as rd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

#set the director for the final file to be written into
os.chdir(r'C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\PlottingImages')

total = 99
perfectMatch = 3/total
pathwayMatch = 96/total
noMatch = 0/total

fig = plt.figure(figsize=(12,1.5))
ax = fig.add_subplot(111)

h1 = ax.barh('', width = perfectMatch, left=0,color='#499562',height=1,align='center') #green, perfect match
h2 = ax.barh('', width = pathwayMatch, left=perfectMatch,color='#FCD238',height=1,align='center') #Yellow, parthway match
h3 = ax.barh('', width = noMatch, left =(perfectMatch + pathwayMatch),color='#CF3B2A',height=1,align='center') #red, no match
             
handles = [h1, h2, h3]
        
ax.axes.get_yaxis().set_visible(False)
ax.xaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
ax.axes.get_yaxis().set_visible(False)

h1 = plt.text( perfectMatch/2, 0., '3', ha="center", va="center", color="white", fontsize=18, fontweight="bold")

plt.text( perfectMatch + pathwayMatch/2, 0., '96', ha="center", va="center", color="white", fontsize=18, fontweight="bold")

plt.text( perfectMatch + pathwayMatch + noMatch/2., 0., '', ha="center", va="center", color="white", fontsize=18, fontweight="bold")

fig.subplots_adjust(bottom = 0.5, right = 0.8)


legend_elements = [Patch(facecolor='#499562', edgecolor='#499562',label='Perfect Match'), Patch(facecolor='#FCD238', edgecolor='#FCD238',label='Pathway Match'), Patch(facecolor='#CF3B2A', edgecolor='#CF3B2A',label='False Positive') ]



box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.3,box.width, box.height * 0.9])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.legend(handles=legend_elements, loc='lower center', ncol = 3, fontsize =18, bbox_to_anchor =(0.49, -2))
plt.rcParams.update({'font.size': 18})

plt.show()

os.chdir(r'C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\PlottingImages')
       
#fig.savefig( 'scale_Phen.png', format='png', dpi=2400)   

#####################################################################
#####################################################################

#write out the structures for the false positives. Use a image program to compile these images with waterfallPlot_LitMatches.py
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import DrawingOptions
import xlsxwriter

#This code can be used to re-write the smiles code as one that will match the others in the network
oldS = 'OC(=O)C1=CC=CC=C1C(O)=O'
newS = Chem.MolToSmiles(Chem.MolFromSmiles(oldS),True) 
print(newS)

###*******************************************************************
DrawingOptions.atomLabelFontSize = 55
DrawingOptions.dotsPerAngstrom = 100
DrawingOptions.bondLineWidth = 9.0

DrawingOptions.padding = 0
DrawingOptions.dblBondOffset = 0.4
DrawingOptions.atomLabelFontFace = 'Times-Bold'

DrawingOptions.elemDict= {0: (0.5, 0.5, 0.5), 1: (0.55, 0.55, 0.55), 7: (0, 0, 1), 8: (0, 0, 0), 9: (0.2, 0.8, 0.8), 15: (1, 0.5, 0), 16: (0.8, 0.8, 0), 17: (0, 0.8, 0), 35: (0.5, 0.3, 0.1), 53: (0.63, 0.12, 0.94)}

os.chdir('C:\ResearchWorkingDirectory\TempImages')

m = Chem.MolFromSmiles('[O-]C(=O)C(=O)c1ccccc1C([O-])=O')
imageSTR = 'Phen1.png'

#comment this out if you don't want the file for the metabolite image to be written out. 
#Chem.Draw.MolToFile(m, imageSTR, size=(700, 700),fitImage=True)

            