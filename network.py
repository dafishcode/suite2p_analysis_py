#=======================================================================
def spacek(coordlist, Fdrop, experiment, mcc):    # Spatial K-means clustering on cell coordinates 
#======================================================================= 
# This function takes the x-y dimensions from identified cells in a numpy array 
# It then performs clustering to pull out spatially contiguous groups of cells  
# mcc = mean cells per cluster
#each plane is clustered separately


#TEST

    import numpy as np
    import os
    from sklearn.cluster import KMeans

    klist = list(range(len(coordlist)))
    for y in range(len(coordlist)):
        print('Clustering fish ' + str(y + 1)+ ' of ' + str(len(coordlist)))
    # Pull out coordinates and loop through each plane
        kvector = []
        cs = np.load(Fdrop + 'Project/' + experiment + os.sep + coordlist[y])
        n_clust  = int(cs.shape[0]/mcc) #how many clusters to make
        kmeans   = KMeans(n_clusters=n_clust, random_state=0).fit(cs[:,:2])  #perform k means on all cells
        kvector =  np.append(kvector, kmeans.labels_) #vector of all label values
        kcoordcs = np.column_stack((cs, kvector)) #array for new coordinates including spatial clusters
        klist[y] = kcoordcs
        np.save(Fdrop + 'Project/' + experiment + os.sep + coordlist[y][:coordlist[y].find('run')+6] + '_' + 'kcoord.npy', kcoordcs)
    return(klist)
    
#=======================================================================
def funck(coordlist, Fdrop, experiment, mcc):    # Spatial K-means clustering on cell coordinates 
#======================================================================= 
    

    
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
        label  = np.load(coordlist[y])[:,3] #cluster labels
    
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
def divconq(kl, cs):   # K-labels, coordinates
#======================================================================= 
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    import fish_net_functions as fn
    from sklearn.cluster import KMeans
    
    kls = np.unique(kl)
    mxdia = []
    for k in kls:
        kid    = np.where(kl == k)[0]
        mxd    = np.max(cs[kid,:] - np.mean(cs[kid,:], axis=0))
        mxdia  = np.append(mxdia, mxd)
        
    toolongs = np.where(mxdia > 100)[0]
    
    if toolongs.shape[0] == 0:     
        return kl
    
    else:  
        nkc = np.max(kls)
        nkl = copy.deepcopy(kl)
        for kcheck in toolongs:
            okc = kcheck
            nkc = nkc + 1          # new counter

            kmembs   = np.where(kl == okc)[0]
            kmeans   = KMeans(n_clusters=2, random_state=0).fit(cs[kmembs,:])
            nkl[kmembs[np.where(kmeans.labels_ == 1)[0]]] = nkc 
        
        nkl = fn.fish_net_divconq(nkl,cs)
        return nkl