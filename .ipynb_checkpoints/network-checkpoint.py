#PROCESS
#------------
#------------

#=======================================================================
def spacek(coordlist, Fdrop, experiment, mcc):    # Spatial K-means clustering on cell coordinates 
#======================================================================= 
# This function takes the x-y dimensions from identified cells in a numpy array 
# It then performs clustering to pull out spatially contiguous groups of cells  
# mcc = mean cells per cluster
#each plane is clustered separately

    import numpy as np
    import os
    from sklearn.cluster import KMeans

    klist = list(range(len(coordlist)))
    for y in range(len(coordlist)):
        print('Clustering fish ' + str(y + 1)+ ' of ' + str(len(coordlist)))
    # Pull out coordinates and loop through each plane 
        kvector = []
        cs = np.load(Fdrop + 'Project/' + experiment + os.sep + coordlist[y])[:,:3]
        spatial_conversion = [.5, .5, 15]
        spacecs = np.multiply(cs, spatial_conversion)
        n_clust  = int(cs.shape[0]/mcc) #how many clusters to make
        kmeans   = KMeans(n_clusters=n_clust, random_state=0).fit(spacecs)  #perform k means on all cells
        kvector =  np.append(kvector, kmeans.labels_) #vector of all label values
        kcoordcs = np.column_stack((cs, kvector)) #array for new coordinates including spatial clusters
        klist[y] = kcoordcs
        np.save(Fdrop + 'Project/' + experiment + os.sep + coordlist[y][:coordlist[y].find('run')+6] + '_' + 'realcoord.npy', kcoordcs)
    return(klist)

#=======================================================================
def funck(Fdrop, experiment, ktrace, kcoord):    # K-means clustering on correlation matrix
#======================================================================= 
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    from sklearn.cluster import KMeans
    import os
    
    #K means cluster correlation matrix
    kcoordlist = list(range(len(ktrace)))
    for y in range(len(ktrace)):
        print('Clustering fish ' + str(y + 1)+ ' of ' + str(len(ktrace)))
        corr = np.corrcoef(np.load(ktrace[y])) 
        coord = np.load(kcoord[y])[:,:3]
        kmeans = KMeans(n_clusters=16, random_state=0).fit(corr)
        klabel = kmeans.labels_
        kcoordnew = np.column_stack((coord,klabel)) #add labels to original kcoord array
        kcoordlist[y] = kcoordnew
        np.save(Fdrop + 'Project/' + experiment + os.sep + kcoord[y][:kcoord[y].find('run')+6] + '_' + 'kcoord.npy', kcoordnew)
    return(kcoordlist)


    
#=======================================================================
def average(Fdrop, experiment, tracelist, coordlist): 
#======================================================================= 
    import numpy as np
    import os
    meantracelist = list(range(len(coordlist)))
    meanloclist = list(range(len(coordlist)))
    for y in range(len(coordlist)):
        print('Calculating fish ' + str(y + 1)+ ' of ' + str(len(coordlist)))
        trace = np.load(tracelist[y]) #trace for each cell
        loc = np.load(coordlist[y])[:,:3] #cell coordinates
        label  = np.load(coordlist[y])[:,np.load(coordlist[y]).shape[1]-1] #cluster labels
    
        labels = np.unique(label) #unique cluster labels
        count     = 0 
    
        for tl in labels: #loop through unique clusters
            cluster      = np.where(label == tl)[0] #all cells in current cluster
            meantrace  = np.mean(trace[cluster,:], axis=0) #average trace for all cells in cluster
            meantracearr     = meantrace if count == 0 else np.vstack((meantracearr, meantrace)) #if first cluster, array first row is meantrace, else vstack each cluster trace 
            
            locmean = np.mean(loc[cluster,:], axis=0) #average location of all cells in cluster
            mloc = locmean if count == 0 else np.vstack((mloc,locmean)) #stack cluster means in vector 

            count  = count + 1
        np.save(Fdrop + 'Project/' + experiment + os.sep + coordlist[y][:coordlist[y].find('run')+6] + '_' + 'ktrace.npy', meantracearr)
        np.save(Fdrop + 'Project/' + experiment + os.sep + coordlist[y][:coordlist[y].find('run')+6] + '_' + 'kcoord.npy', mloc)
        meantracelist[y] = meantracearr
        meanloclist[y] = mloc
    return(meantracearr, mloc)



#======================================================================= 
def divconq(kcoord, i, Fdrop, experiment, kcoordinput, kread):   # K-labels, coordinates
#======================================================================= 
    import copy
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    from sklearn.cluster import KMeans
    
    cluscoord = kread[:,:3] #spatial cluster coordinates
    kl = kread[:,3] #functional labels of spatial clusters
    kls = np.unique(kl) #unique functional label
        
    mxdia = [] #max distance of s-cluster from its assigned f-cluster for all clusters
    for k in kls: #loop through unique label values
        kid    = np.where(kl == k)[0] #indeces of current label values
        mxd    = np.max(cluscoord[kid,:] - np.mean(cluscoord[kid,:], axis=0)) #max distance from mean of cluster
        mxdia  = np.append(mxdia, mxd) 
        
    toolongs = np.where(mxdia > 100)[0] #returns the value for clusters which are outside distance threshold

    if toolongs.shape[0] == 0:   #if no clusters are above this distance - keep functional labels as before   
        return kl
    else:  #if some cluster are above max distance
        nkc = np.max(kls)     #max label value
        nkl = copy.deepcopy(kl) #copy of kl
        for kcheck in toolongs: #loop through each cluster > 100 in length
            okc = kcheck #current 
            nkc = nkc + 1          #iterate over cluster number max value
            kmembs   = np.where(kl == okc)[0] # where functional label == current label loop
            kmeans   = KMeans(n_clusters=2, random_state=0).fit(cluscoord[kmembs,:]) #repeat kmeans - and split into 2 smaller cluster
            nkl[kmembs[np.where(kmeans.labels_ == 1)[0]]] = nkc #assign new cluster values
            kcoordnew = np.column_stack((cluscoord,nkl))
        nkl = divconq(kcoord, i, Fdrop, experiment, kcoordinput, kcoordnew) #iterate over all remaining clusters until no cluster distance > 100
    print('Clustered fish ' + str(i + 1)+ ' of ' + str(len(kcoord)))
    return(nkl)


#===============================================================================
def watts_strogatz(edge_density, p, Nnodes, dist, savepath): 
#==============================================================================

    """Generate random small world graph with specific Edge density. The Watts-Strogatz model has (i) a small average shortest path length, and (ii) a large clustering coefficient. The algorithm works by assigning a pre-defined number of connections between k-nearest neighbours - it then loops through each node and according to some uniform probability re-assigns its edges from its connected k-nearest neighbours and a random unconnected node. 

    edge_density = number of k_nearest neighbours each node is connected to
    p = probability of an edge being randomly re-assigned
    Nnodes = number of nodes
    dist = distance matrix between all nodes in network
    """


    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import networkx as nx
    from sklearn.metrics.pairwise import euclidean_distances
    import copy
    import random
    
    k_neighbours  = int(Nnodes * edge_density)
    A  = np.zeros(dist.shape)

    if k_neighbours == 0: 
        return A
    
    # Loop through rows of distance matrix to find k_neighbours
    #-----------------------------------------------------------------------------
    for row in range(dist.shape[0]):
        neighbours = dist[row,].argsort()[:k_neighbours+1][::-1] #find neighbours 
        A[row,neighbours[:-1]] = 1 #make all edges of neighbours connected in network
        A[neighbours[:-1],row] = 1

    # Rewire connections with certain probability
    #-----------------------------------------------------------------------------
    [rows, cols]    = np.where(np.triu(A) == 1) 
    probs           = np.random.uniform(size = rows.shape[0]) #Generate random values for each connection 
    edges_to_change = np.where(probs <= p)[0] #see which values are randomly changed
    old_A  = copy.deepcopy(A) #create copy of A


    for e in range(edges_to_change.shape[0]): #Loop through edges to change
        this_edge = edges_to_change[e]
        A[rows[this_edge], cols[this_edge]] = 0         # switch off old edge
        A[cols[this_edge], rows[this_edge]] = 0

        where_0 = np.where(A[rows[this_edge]] == 0)[0] #find possible connections to reassign to
        new_edge = random.choice(where_0[np.where(where_0 !=rows[this_edge])[0]]) #randomly choose one - ignoring any connections on the diagonal 
        #Assign connection
        A[rows[this_edge], new_edge] = 1        # switch off old edge
        A[new_edge, rows[this_edge]] = 1
    np.save(savepath + 'netparam-k' + str(edge_density) + '-p' + str(p) +  '.npy', old_A)
    return A, old_A

#ANALYSIS
#------------
#------------