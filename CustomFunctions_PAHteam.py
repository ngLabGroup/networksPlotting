
"""
Created on Mon Aug 27 12:05:58 2018

@author: Trevor 
"""
#This contains support functions for the other scripts

def FindSourceSink(G):
    #**********finds all the source nodes and sink nodes****************
    SourceN = []
    SinkN = []
    
    #we know what the source and sink nodes are, but let's write some code to easily find them
    bb = dict(G.out_degree())
    cc = dict(G.in_degree())
    
    for key in bb:
        if bb[key]==0:
#            print ('there is nothing out at ', key)
            SinkN.append(key)
        elif cc[key]==0:
#            print ('there is nothing in at,', key)
            SourceN.append(key)
    #**********finds all the source nodes and sink nodes****************        
    return([SourceN, SinkN])     
    

#
#    '''If there is a cycle that is reachable from root, then result will not be a hierarchy.
#
#       G: the graph
#       root: the root node of current branch
#       width: horizontal space allocated for this branch - avoids overlap with other branches
#       vert_gap: gap between levels of hierarchy
#       vert_loc: vertical location of root
#       xcenter: horizontal location of root            
#    '''
def hierarchy_posMS(G, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5 ):
#    from FindSourceSink import FindSourceSink
    
    
    def h_recur(G, root, width, xcenter, vert_gap = 0.2, vert_loc = 0,  
                      pos = None, parent = None, parsed = [] ):
                 
            if(root not in parsed):
                parsed.append(root)
                if pos == None:
                    pos = {root:(xcenter,vert_loc)}
                else:
                    pos[root] = (xcenter, vert_loc)
                neighbors = list(G[root])
    
                if len(neighbors)!=0:
                    dx = width/len(neighbors)
                    nextx = xcenter - width/2 - dx/2
                    for neighbor in neighbors:
                        nextx += dx
                        #recurse the function except call with "neighbor"
                        pos = h_recur(G,neighbor, width = dx, xcenter=nextx,vert_gap =                  vert_gap, vert_loc = vert_loc-vert_gap,  pos=pos, parent = root, parsed = parsed)
            
            return pos
       
    #Starting the actual code here, everything above is function definitions
    
##temporary code for testing
#    G=nx.DiGraph()
#    G.add_edges_from([(1,2), (1,3), (1,4), (2,5), (2,6), (2,7), (3,8), (3,9), (4,10),(5,11), (5,12), (6,13), (15,2), (15,17), (17,18), (15,10), (17,20), (17,21), (30,20)])
#    
    #locate the roots and sinks based on degrees in/out
    [All_roots, All_sinks] = FindSourceSink(G)
#    print('the roots are', All_roots)
    widthIn = 1/len(All_roots)
    xcenter = widthIn/2
    posfinal = {}
    
    #wrap the recursing code in a loop
    #**************
    for currentRoot in All_roots:
        #each h_recur iterates through a tree
        pos = h_recur(G, currentRoot, width=widthIn, xcenter = xcenter, vert_gap = 0.2, vert_loc = 0)
              
        #append pos onto posfinal
        posfinal = {**posfinal, **pos}
        #clean pos for the next iteration
        pos = {}
        #update xcenter for the next iteration
        xcenter = xcenter + widthIn
    #**************
    return posfinal

    
#Alternative hierarchy plotting
def hierarchy_posLanes(G, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5 ):
    #first thing is to go through the network and assign the "y" values
    from FindSourceSink import FindSourceSink
    [SourceN, SinkN] = FindSourceSink(G)

    #**************
    return posfinal


def getIncomingNodes(G, node, SourceN):
#will return all of the nodes that lie on a pathway between sourceN and node
#******************************def getIncomingNodes(G, node, SourceN):
#will return all of the nodes that lie on a pathway between sourceN and node
#******************************************
#initialize the allPaths variable if on the first call to this function
#node = 'Oc1ccccc1C([O-])=O'
    
    currentS = []
    allSources = []

    inEdges = list(G.in_edges(node))
    for i in inEdges:
        currentS.append(i[0])
    allSources = currentS        
        
    while not(set(currentS).issubset(set(SourceN))):
    #    allPaths = newPaths
        #build a array of paths
        newPaths = [] 
    
        #next round of sources to check and see if we're done or not. 
    #    currentS = []
        nextS = []
        #need to change this
        for c in currentS:
            #p[0] is the most recent node in the list
            inEdges = list(G.in_edges(c))
    
            
            #These are the new edges at the top of p
            for x in inEdges:
                #This is appending even if the path is already there
                #let's try just keeping track of the nodes, not all the edges also. 
                
    #            newPaths.append([x[0]]+p)
                nextS.append(x[0])
                
                
        currentS = list(set(nextS))       
        allSources = list(set(allSources + nextS))
#        print(str(len(currentS) ) + ' and overall length is ' + str(len(set(allSources))))
    
    return allSources
    #******************************************************