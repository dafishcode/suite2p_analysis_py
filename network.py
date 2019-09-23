#=======================================================================
def kmeans(coordlist, Fdrop, experiment, mcc):    # K-means clustering on cell coordinates 
#======================================================================= 
# This function takes the x-y dimensions from identified cells in a numpy array 
# It then performs clustering to pull out spatially contiguous groups of cells  
# mcc = mean cells per cluster

    import numpy as np
    import os
    from sklearn.cluster import KMeans
    
    # Pull out coordinates and loop through each plane
    cs = np.load(Fdrop + 'Project/' + experiment + os.sep + coordlist[0])
    
    for pl in np.unique(cs[:,2]):  # loop through each plane, third dimension contains the plane index 
        print('Now working on plane ' + str(pl) + ' of ' + str(np.max(np.unique(cs[:,2]))))
        id       = np.where(cs[:,2] == pl)[0] #locations for cells of current plane
        n_clust  = int(np.max(id.shape) / mcc)
        kmeans   = KMeans(n_clusters=n_clust, random_state=0).fit(cs[id,0:2])
        
                 
                 
        Fish["Planes"][int(pl)].update({"Klabel":kmeans.labels_, "Index":id})
   
    # Assemble everthing into single k-label list
    #---------------------------------------------------------------------
    ks = np.empty(0)
    i  = 0

    for p in Fish["Planes"]:
        mx = 0 if i == 0 else np.max(ks)
        ks = np.append(ks, p["Klabel"] + mx)
        i = i + 1
    
    Fish.update({"Klabel":ks})
    return Fish
    