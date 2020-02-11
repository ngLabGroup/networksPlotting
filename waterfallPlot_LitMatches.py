# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 12:38:20 2019

@author: twsle
"""
get_ipython().magic('reset -sf')


#************************************************
#Waterfall Plots with literature Matchs
#This script creates the waterfall style plot with the compounds matched 
#by empirical literature plotted as images rather than yellow dots

#************************************************

import os
import gc
gc.collect()

os.chdir(r'C:\ResearchWorkingDirectory\gitRepositories\networksPlotting')
from CustomFunctions_PAHteam import FindSourceSink
from CustomFunctions_PAHteam import hierarchy_posMS

import networkx as nx
import pandas as pd
import numpy as np
import xlsxwriter

from rdkit import DataStructs
from rdkit import Chem

from rdkit.Chem.Draw import DrawingOptions
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
from PIL import ImageFile

from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

DrawingOptions.atomLabelFontSize = 80
DrawingOptions.dotsPerAngstrom = 900
DrawingOptions.bondLineWidth = 11.0
DrawingOptions.padding = 0
DrawingOptions.dblBondOffset = 0.4
DrawingOptions.atomLabelFontFace = 'Times-Bold'

#8 is the oxygens. This sets all oxygens to black. 
DrawingOptions.elemDict= {0: (0.5, 0.5, 0.5), 1: (0.55, 0.55, 0.55), 7: (0, 0, 1), 8: (0, 0, 0), 9: (0.2, 0.8, 0.8), 15: (1, 0.5, 0), 16: (0.8, 0.8, 0), 17: (0, 0.8, 0), 35: (0.5, 0.3, 0.1), 53: (0.63, 0.12, 0.94)}

#one could easily write a short segment of code to select the high throughput compounds based on a threshold. However, since it's a small number and it's usually desireable to be able to easily view this dataset I find it most helpful to have them in another excel file. 
df = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\LiteratureMatches_HighThroughput.xlsx", sheet_name = 'HighThroughput')

#select the desired node list.  
#****************************
nodeList = list(df['ExampleSmiles'])
#nodeList = list(df['AnthraceneSmiles'])


#remove the NaNs
nodeList = [nodeList for nodeList in nodeList if str(nodeList) != 'nan']


#Option to remove any high throughput compounds that are not part of a pathway match with a literature node
#****************************************
##phenanthrene
#nodeList.remove('[O-]C(=O)C(=O)c1ccccc1C([O-])=O')
#****************************************

df_raw = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\RawExampleData.xlsx", sheet_name = 'Sheet1')

#need to make sure we use the weights here
G_raw = nx.from_pandas_edgelist(df_raw, source ='From', target = 'To', edge_attr=True, create_using=nx.DiGraph())

newDF = df_raw.iloc[0:0]
allNodes = []
outerNodes = []
innerNodes = []

for n in nodeList:
    
    outEdges = list(G_raw.out_edges(n))
    
    for o in outEdges: 
        
        if (o[0] in nodeList) &  (o[1] in nodeList):
            
            #need to figure out a way to get just the applicable edges, not the junk
            tempc = df_raw.loc[(df_raw['From'] == o[0]) & (df_raw['To'] == o[1])]
            newDF = newDF.append(tempc)
            tempc = []
                        
        else:
            pass
#            print("not in the list ", n)
  
#G_new is a network build from perfect matches and pathway matches        
G_new = nx.from_pandas_edgelist(newDF, source ='From', target = 'To', edge_attr=True, create_using=nx.DiGraph())

nodes = G_new.nodes()

os.chdir('C:\ResearchWorkingDirectory\gitRepositories')

pos = hierarchy_posMS(G_new)

#comment this out if you don't want the pos to reset back to what the code generates. this lets you move compounds manually sothat you can get a better looking plot. 
posOut = pd.DataFrame.from_dict(pos, orient = 'index')
posOut.to_excel('ExamplePositions.xlsx', sheet_name='Sheet1', engine='xlsxwriter' )

#This file manually moves the compounds around so that you can see the compounds. 
##phenanthrene
posIn = pd.read_excel(r"C:\ResearchWorkingDirectory\gitRepositories\ExamplePositions.xlsx", sheet_name = 'Sheet1')

posIn = posIn[[0,1]]
pos = {}
for pIn in posIn.iterrows():
    pos[pIn[0]] = list(pIn[1])

#This is the literature data, nodes in this list match literature nodes
knownNodesDF = pd.read_excel(r"C:\ResearchWorkingDirectory\DataReferenceFiles\PaperData\SimiliarityCompounds.xlsx", sheet_name = 'Phenanthrene')
listNodes = list(knownNodesDF['Perfect_Match'])
perfectMatch = [x for x in listNodes if str(x) != 'nan']

#initialized the variables we need for the next part. 
colorList = []
nodes1 = []
nodes2 = []
nodes3 = []

for n in nodes:
    if n in perfectMatch:
        #perfect match - green color if you're not plotting these as images
        colorList.append('#499562')
        nodes1.append(n)
    else:
        #pathway match - yellow color
        colorList.append('#FCD238')
        nodes2.append(n)
        #Also build the labels for these at the same time
        
        
#get the edges for the perfect match and the pathway matches
fullEdgelist = list(G_new.edges())
litMatches_edgelist = []
pathMatches_edgelist = []
for t in fullEdgelist:
    #select the edges that come to or from a 
    pass
    if (t[1] in perfectMatch ):
        litMatches_edgelist.append(t)
    else:
        pathMatches_edgelist.append(t)

##draws both nodes and edges, the marker is supposed to cover up enough of the edge that you can see the image structure, only use the literature matches because we have the big nodes here
nc = nx.draw_networkx(G_new, pos, marker = "o", node_color = 'white', with_labels=False,  node_size=3400, nodelist=nodes1, edgelist = litMatches_edgelist, zorder = 1) # vmin = 0, vmax = 6)


#draw the pathway matches, use the pathways match edgelist
ec = nx.draw_networkx_edges(G_new, pos, edgelist = pathMatches_edgelist, alpha=0.8, width = 1, zorder = 3)
  
##Pathway Match
nc = nx.draw_networkx_nodes(G_new, pos, marker = "o", node_color = '#FCD238', edge_color = None, nodelist=nodes2, labels =list(range(1, len(nodes2))) ,  node_size=200, zorder = 5) # vmin = 0, vmax = 6)                                                   

                           
iterator = 1
labelDF = pd.DataFrame(columns=['SMILES', 'x', 'y'])

for n in nodes2:
    tempSlice = pd.DataFrame([{'SMILES': n, 'x': pos[n][0], 'y':pos[n][1]}]) 
    labelDF = labelDF.append(tempSlice) 
    iterator +=1
    
#****   sort the labels top to bottom so that the labels descend top to bottom  *****
labelDF.sort_values(by=['y'], inplace = True, ascending=False,)    
  
uniqueY = list(set(labelDF['y']))
uniqueY.sort(reverse = True)

#initialize labelDF2
labelDF2 = pd.DataFrame(columns=['SMILES', 'x', 'y'])

#Rebuild the labelDF, sorting left to right within each layer as you go 
for u in uniqueY:
    
    tempSlice = labelDF.loc[labelDF['y'] == u]
    tempSlice.sort_values(by=['x'], inplace = True, ascending=True,)
    labelDF2 = labelDF2.append(tempSlice)
      
labelDF2.reset_index(inplace=True, drop = True)  
labelDF2.drop(columns=['x', 'y'], inplace = True)
labelDF2['Labels'] = labelDF2.index

labelDF2.set_index('SMILES', drop = True, inplace = True)

labelDict = {}
for l in labelDF2.iterrows():
    labelDict[l[0]] =     l[1]['Labels']+1

#labelDict is the final product needed to use for labeling the nodes2 section
nx.draw_networkx_labels(G_new,pos, nodelist=nodes2,labels =labelDict,font_size=12)                            
                            
plt.show()
axes = plt.gca()

os.chdir('C:\ResearchWorkingDirectory\gitRepositories\PlottingImages')

#highlight specific nodes if desired
#otherSinks = ['Oc1cc(O)c(O)c(O)c1O']
otherSinks = []

for p in pos:
    m = Chem.MolFromSmiles(p)

    imageSTR = 'image' + str(iterator) + '.png'
    Chem.Draw.MolToFile(m, imageSTR, size=(700, 700),fitImage=True, BoldText = True, fill = None)
        
    img = Image.open('C:\ResearchWorkingDirectory\gitRepositories\PlottingImages\\' + imageSTR)
    img = img.convert("RGBA")
    datas = img.getdata()
    
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    
    img.putdata(newData)
    img.save('C:\ResearchWorkingDirectory\gitRepositories\PlottingImages\\' + imageSTR, "PNG")
     
    imagebox = OffsetImage(img, zoom=.125, zorder = 1)
    
    tempCoords = list(pos[p])

    if p in perfectMatch:        
        ab = AnnotationBbox(imagebox, tempCoords, xybox= (0,-0),xycoords='data',bboxprops = dict(boxstyle=("Square, pad=0.1"), color = 'None',facecolor='None', linewidth=.1), boxcoords="offset points")
        ab.set_zorder(10)
        print('found it')
        axes.add_artist(ab)
    elif p in otherSinks:
        
        ab = AnnotationBbox(imagebox, tempCoords, xybox= (0,-0),xycoords='data',bboxprops = dict(boxstyle=("Square, pad=0.1"), color='None', linewidth=.1), boxcoords="offset points")
       
        ab.set_zorder(10)
        print('matched ',p)
#        axes.add_artist(ab)
    else:
        pass

    iterator += 1
        
#*******************************
#Generate a high quality image

fig = plt.gcf()
cur_axes = plt.gca()
cur_axes.axes.get_xaxis().set_visible(False)
cur_axes.axes.get_yaxis().set_visible(False)
cur_axes.axis('off')

#create an image so that we can crop out the white space and label it if desired
fig.set_size_inches(12, 6)
fig.savefig('myimage1.png', format='png', dpi=2400)
Image.MAX_IMAGE_PIXELS = None
image=Image.open('myimage1.png')
image.load()

image_data = np.asarray(image)
#image_data_bw = image_data.max(axis=2)
#
#
#non_empty_columns = np.where(image_data_bw.max(axis=0)>0)[0]
#non_empty_rows = np.where(image_data_bw.max(axis=1)>0)[0]
#cropBox = (min(non_empty_rows), max(non_empty_rows), min(non_empty_columns), max(non_empty_columns))

#manual crop to waste a little less white space
image_data_new = image_data[1000:13000,5000:25000 , :]

##Option to add some blank space and text on top
#**************************************
#blankSpace = np.full((1000,20000,4), 0, dtype = 'uint8')
#image_data_final = np.vstack( (blankSpace,image_data_new))

new_image = Image.fromarray(image_data_new)

#option to label the plot
#font = ImageFont.truetype("arial.ttf", 550)
#ImageDraw.Draw(new_image).text((100, 200),  'Anthracene High Throughput Predicted Metabolites' , (0, 0, 0) ,font = font)

#new_image.save('Phen_Labels_final.png')
#new_image.close()



