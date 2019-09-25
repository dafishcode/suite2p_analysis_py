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
        kmeans = KMeans(n_clusters=100, random_state=0).fit(corr)
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



#ANALYSIS
#------------
#------------