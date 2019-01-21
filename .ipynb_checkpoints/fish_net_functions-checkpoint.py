#=======================================================================
def fish_net_load(Fdata): # Load imaging datasets to custom object
#=======================================================================    
# This function looks in the Fdata folder for specific files from the suite2p output
# 1) A list of numpy arrays with the plane information - of the form 'plane0... .npy'
# 2) A single numpy file containing the cell coordinates - called 'com_signal... .npy'
# 3) A single numpy file containing a ncells x ntimepoints array of data - called 'com_signal... .npy'

    import os 
    import re
    import numpy as np

    dirlist = os.listdir(Fdata)

    # Find planes of suite2p output
    #---------------------------------------------------------------------
    r       = re.compile('^plane[0-9].*')
    planelist = list(filter(r.match, dirlist))
    planelist.sort()

    P = []
    for p in planelist:
        path = Fdata + os.sep + p
        p_id = int(p[5:p.find('_')])
        P.append({"path":path, "plane_id":p_id})

    print('Found ' + str(len(planelist)) + ' planes') 

    #     P[i]["id"]   = p

    # Find coordinates of suite2p output
    #---------------------------------------------------------------------
    r       = re.compile('^com_coord')
    coord   = list(filter(r.match, dirlist))
    coord   = np.load(Fdata + os.sep + coord[0])


    # Find combined signal trace
    #---------------------------------------------------------------------
    r       = re.compile('^com_signal')
    signl   = list(filter(r.match, dirlist))
    signl   = np.load(Fdata + os.sep + signl[0])

    Fish = {"Planes":P, "Data":signl, "Coordinates":coord}

    return Fish

#========================================================================
def fish_net_regrout(dat, fun = 'lin'):   # Remove slow drift components 
#========================================================================
# This function removes some baseline / slow drift in the fluorescence of the data 
# This can either be done by removing a single linear regression component ('lin')
# or high pass filtering the data so that slower drifts are removed ('filt')

    import numpy as np
    from scipy import stats, signal
    
    alld = np.zeros(dat.shape)
    if fun == 'lin':
        
        for i in range(dat.shape[0]):
            d = dat[i,:]
            slope, intercept, a,b,c = stats.linregress(np.linspace(0,len(d)-1,len(d)),d)  
            d = d - (slope*d + intercept)   
            
            alld[i,:]  = d
            
    elif fun == 'filt':
        
        # Filter specs
        #----------------------------------------------------------------------------------------------
        fc  = .001           # cutoff frequency
        fs  = 2.7            # sampling frequency
        nfc = fc / (fs / 2)  # normalised cutoff 

        # Generate and apply filter
        #----------------------------------------------------------------------------------------------
        b, a   = signal.butter(5, nfc, 'high')
        alld   = signal.filtfilt(b,a, dat, axis=1,method='gust')
        
    
    return alld


#=======================================================================
def fish_net_spacek(Fish, mcc = 10):    # K-means clustering on cell coordinates 
#======================================================================= 
# This function takes the x-y dimensions from identified cells in a numpy array 
# It then performs clustering to pull out spatially contiguous groups of cells
# (This is borrowed from Rick Betzel's zebrafish paper: 
#  https://www.biorxiv.org/content/early/2018/12/15/496414 )    
# 
# mcc = mean cells per cluster
# returns updated Fish object

    import numpy as np
    from sklearn.cluster import KMeans
    
    # Pull out coordinates and loop through each plane
    cs = Fish["Coordinates"]
    
    for pl in np.unique(cs[:,2]):  # third dimension contains the plane index 
        print('Now working on plane ' + str(pl) + ' of ' + str(np.max(np.unique(cs[:,2]))))
        id       = np.where(cs[:,2] == pl)[0]
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
    

    
#=======================================================================
def fish_net_average(d, l, loc): # d = data matrix, l = labels (numeric)
#======================================================================= 
    import numpy as np
    
    labels = np.unique(l)
    ct     = 0
    ml     = []
    
    for tl in labels:
        i      = np.where(l == tl)[0]
        nmean  = np.mean(d[i,:], axis=0) 
        md     = nmean if ct == 0 else np.vstack((md, nmean))
        ml     = np.append(ml,tl)
        
        locmean = np.mean(loc[i,:], axis=0)
        mloc = locmean if ct == 0 else np.vstack((mloc,locmean))

        ct  = ct + 1

    return md, ml, mloc



#=======================================================================
def fish_net_nneigh(cs, rng = 6000, dim = [.8, .8, 15], cnt=5): # xyz (or xy) coordinates of nodes
#======================================================================= 
    import numpy as np
    
    # Set up nearest neighbour graph
    #---------------------------------------------------------------------------
    mcs  = np.multiply(cs, dim)     # metrically scaled coordinates (in microns)

    # Initialise full distance matrix and nearest neighbour graph (binary) matrix
    #---------------------------------------------------------------------------
    nnb  = np.zeros((cs.shape[0],cs.shape[0]))

    # Loop through all cells to fill in distances
    #---------------------------------------------------------------------------
    for r in range(cs.shape[0]):
        dis = np.ones((10,cs.shape[0]))*10000

        if r % round((10*cs.shape[0]/100)) == 0: print("Doing row " + str(r) + " of " + str(cs.shape[0])) 
        for rr in range(max([r-int(rng/2),0]), min([r+int(rng/2),dis.shape[1]])):  # moving window around r

            if r == rr: dis[0,rr] = 10000 
            else:       dis[0,rr] = np.linalg.norm(mcs[r,:]-mcs[rr,:])

        mini = np.where(dis[0,:] < np.nanpercentile(dis[0,:],cnt))[0]
        nnb[r,mini] = 1

    print('Done')
    
    return nnb

#======================================================================= 
def fish_net_peaks(dat, cnt = 95, typ = 'std', stdlim = 3):
#======================================================================= 
    import numpy as np
    from scipy import stats
    from scipy import signal
    
    # Find activity peaks
    #---------------------------------------------------------------------------
    pks = np.zeros(dat.shape)
    for i in range(dat.shape[0]):
        d = dat[i,:]                                                            
        if typ == 'peaks':
            p, prop = signal.find_peaks(d,threshold=np.percentile(d,cnt))
        elif typ == 'std':
            sem = np.std(d)
            p   = np.where(d > stdlim*sem)[0]
        else: print('Don''t know what type of binarisation to use')
        
        for pi in p: pks[i,p] = 1
            
        
    
    return pks

#======================================================================= 
def fish_net_avalanche(pks, nnb): 
#======================================================================= 
# pkg - returns peaks grouped into contiguous clusters (coded by numbers)
# avsz - returns the size of all avalanches (in integers)

    import numpy as np
    
    pkg    = np.zeros(pks.shape)                       # peak groups
    act_t  = np.where(np.sum(pks, axis=0) > 3)[0]      # Time points with at least 3 cells active

    i = 0 
    for t in act_t:
        
        tlen = act_t.shape[0]
        if i % round(10*tlen/100) == 0: print("Doing time point " + str(i) + " of " + str(tlen)) 
        i   = i + 1
        
        gr  = 1
        cid = np.where(pks[:,t] > 0)[0]
        
        for c in cid:
            
            # If currently unlabelled, label with gr + 1
            #-------------------------------------------------------------------
            if pkg[c,t] == 0:        
                gr = gr + 1
                pkg[c,t] = gr

            # Find all neighbours
            #-------------------------------------------------------------------
            nb   = np.where(nnb[c,:] > 0)[0]
            tgr  = np.intersect1d(cid,nb)


            # Label those that are active the same as 
            #-------------------------------------------------------------------
            pkg[tgr,t] = pkg[c,t]
            
    # For each time point count number of coactivations
    #----------------------------------------------------------------------------
    print('Now calculating avalanche size')
    avsz = np.array([])
    for t in range(pkg.shape[1]):
        comms = np.unique(pkg[:,t][pkg[:,t] > 0])
        for c in comms:
            avsz = np.append(avsz, pkg[:,t][pkg[:,t] == c].shape[0])    

    return pkg, avsz


#======================================================================= 
def fish_net_divconq(kl, cs):   # K-labels, coordinates
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