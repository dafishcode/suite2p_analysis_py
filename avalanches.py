import admin_functions as adfn
import IS as isfn


#PROCESS
#------------
#------------
#=======================================================================
def neighbour_new(coord, n_neigh, dim): # Select which fish data to visualise
#=======================================================================
    import numpy as np
    import os
    
    #Loop through all fish
    #----------------------
        
        # Set up nearest neighbour graph
        #---------------------------------------------------------------------------
    mcs  = np.multiply(coord, dim)     # metrically scaled coordinates (in microns)
        
        # Initialise full distance matrix and nearest neighbour graph (binary) matrix
        #nearest neigh binary matrix of celln by celln storing 
        #distance of each cell to every other cell
        #---------------------------------------------------------------------------
    nnb  = np.zeros((coord.shape[0],coord.shape[0]))  
        
    for r in range(coord.shape[0]):
        distance = np.zeros(coord.shape[0])
        if r % round((10*coord.shape[0]/100)) == 0: 
            print("Doing row " + str(r) + " of " + str(coord.shape[0]))
        
        for x in range(coord.shape[0]):
            if x == r: 
                distance[x] = 100000
            else:
                distance[x] = np.linalg.norm(mcs[r]-mcs[x]) 
                
        index = np.argsort(distance)[:n_neigh]
        nnb[r,index] = 1 #binary value defining whether in range or not 
    return(nnb)

#=======================================================================
def neighbour_r(coord, cnt, rng, dim): # Select which fish data to visualise
#=======================================================================
    import numpy as np
    import os
    
    #Loop through all fish
    #----------------------
        
        # Set up nearest neighbour graph
        #---------------------------------------------------------------------------
    mcs  = np.multiply(coord, dim)     # metrically scaled coordinates (in microns)
        
        # Initialise full distance matrix and nearest neighbour graph (binary) matrix
        #nearest neigh binary matrix of celln by celln storing 
        #distance of each cell to every other cell
        #---------------------------------------------------------------------------
    nnb  = np.zeros((coord.shape[0],coord.shape[0]))  
        
        # Loop through all cells to fill in distances
        #distance = matrix of celln x planen *10000 so default value is v large, 
        #outside of typical range and then will fill with distances for connected cells
        #---------------------------------------------------------------------------
    for r in range(coord.shape[0]):
        distance = np.ones((10,coord.shape[0]))*10000
        if r % round((10*coord.shape[0]/100)) == 0: 
            print("Doing row " + str(r) + " of " + str(coord.shape[0]))
            
            # moving window around r of size 3000 cells either side 
            # for each value of cell(r), each rr value (cell that is within range of cellr) 
            # a distance is calculated from cell to rrcell from their metrically scaled positions in space
            #------------------------------------------------------------------------------------
        #CHANGE TO SELECTION OF N CLOSEST NEURONS - SAME N ACROSS FISH
        for rr in range(max([r-int(rng/2),0]), min([r+int(rng/2),distance.shape[1]])):  
            if r == rr: distance[0,rr] = 10000  #cannot connect to itself - set to 10000 ie value too large to be in range 
            else:       distance[0,rr] = np.linalg.norm(mcs[r,:]-mcs[rr,:]) 
            
            #calculate binary matrix of all cells that are in range
            #--------------------------------------------------------------
        mini = np.where(distance[0,:] < np.nanpercentile(distance[0,:],cnt))[0] #THIS IS THROWN OFF BY ALL THE 10000S - JUST CONVERT TO SELECTING NUMBER OF NEURONS
        nnb[r,mini] = 1 #binary value defining whether in range or not 
    return(nnb)


#=======================================================================
def neighbour(cnt, savepath, experiment, rng, dim, name): # Select which fish data to visualise
#=======================================================================
    import numpy as np
    import os
    #Loop through all fish
    #----------------------
    coord = np.load(name)[:,:3]
        
        # Set up nearest neighbour graph
        #---------------------------------------------------------------------------
    mcs  = np.multiply(coord, dim)     # metrically scaled coordinates (in microns)
        
        # Initialise full distance matrix and nearest neighbour graph (binary) matrix
        #nearest neigh binary matrix of celln by celln storing 
        #distance of each cell to every other cell
        #---------------------------------------------------------------------------
    nnb  = np.zeros((coord.shape[0],coord.shape[0]))  
        
        # Loop through all cells to fill in distances
        #distance = matrix of celln x planen *10000 so default value is v large, 
        #outside of typical range and then will fill with distances for connected cells
        #---------------------------------------------------------------------------
    for r in range(coord.shape[0]):
        distance = np.ones((10,coord.shape[0]))*10000
        if r % round((10*coord.shape[0]/100)) == 0: 
            print("Doing row " + str(r) + " of " + str(coord.shape[0]) + " for " + name[:name.find('sess')-9] + '_' + name[name.find('dpf')+4:name.find('run')-1])
            
            # moving window around r of size 3000 cells either side 
            # for each value of cell(r), each rr value (cell that is within range of cellr) 
            # a distance is calculated from cell to rrcell from their metrically scaled positions in space
            #------------------------------------------------------------------------------------
        for rr in range(max([r-int(rng/2),0]), min([r+int(rng/2),distance.shape[1]])):  
            if r == rr: distance[0,rr] = 10000  #set to 10000 ie value too large to be in range
            else:       distance[0,rr] = np.linalg.norm(mcs[r,:]-mcs[rr,:]) 
            
            #calculate binary matrix of all cells that are in range
            #--------------------------------------------------------------
        mini = np.where(distance[0,:] < np.nanpercentile(distance[0,:],cnt))[0]
        nnb[r,mini] = 1 #binary value defining whether in range or not 
    np.save(savepath + 'Project/' + experiment + os.sep + name[:name.find('run')+6] + '_' +    str(cnt) + 'nnb.npy', nnb)
    return(nnb)


#=======================================================================
def corrdis_bin(corr, dist, bins):
#=======================================================================
    import numpy as np
    if corr.shape[0] != dist.shape[0]:
        print('Correlation and Distance matrices have unequal cell numbers')
        return()
    
    # Take upper triangular of matrix and flatten into vector
    corr = np.triu(corr, k=0) 
    dist = np.triu(dist, k=0)
    corr_v = corr.flatten()
    dist_v = dist.flatten()

    # Convert all negative correlations to 0
    corr_v = [0 if o < 0 else o for o in corr_v]
    corr_v = np.array(corr_v)
    dist_v[np.where(corr_v == 0)] = 0

    # Order by distances
    unq = np.unique(dist_v)
    dist_vs = np.sort(dist_v)
    corr_vs = np.array([x for _,x in sorted(zip(dist_v,corr_v))])
    window = adfn.window(np.int((unq.shape[0]/bins)), unq.shape[0])[0] #define a sliding window

    #Loop through each bin and calculate average distance/correlation
    count, bincount=0,0
    dist_bins, corr_bins = np.zeros(np.int(unq.shape[0]/window)),np.zeros(np.int(unq.shape[0]/window))
    for i in range(np.int(unq.shape[0]/window)):
        if i == np.int(unq.shape[0]/window)-1:
            break
        start = count
        stop = count+window
        start_in = np.where(dist_vs == unq[start])[0][0] 
        stop_in = np.where(dist_vs == unq[stop])[0][0] - 1

        sumd_c = np.sum(corr_vs[start_in:stop_in])
        div_c = len(np.where(corr_vs[start_in:stop_in] !=0)[0])
        corr_bins[bincount] = sumd_c/div_c

        sumd_d = np.sum(dist_vs[start_in:stop_in])
        dist_bins[bincount] = sumd_d/div_c

        bincount+=1
        count+=window
    return(np.vstack((dist_bins, corr_bins)))
    
    
    

#ANALYSIS
#------------
#------------
#=======================================================================
def avalanche_r(nnb, bind): 
#=======================================================================
    import numpy as np
    import os
    import itertools


#Calculate avalanche size + duration
#-----------------------------------
    binarray, oldav, firstav, realav, timemachine, convertav, fill, time = [],[],[],[],[],[],[],[]
    
    #LOOP THROUGH EACH FISH
    #---------------------------------
    #---------------------------------
    binarray, nnbarray, pkg = bind,nnb, np.zeros(bind.shape)
    i, marker, avcount = 0,0,0
        
    #LOOP THROUGH EACH TIME POINT
    #------------------------------
    #------------------------------
    for t in range(binarray.shape[1]-1): #loop through all time points
        #if i% round(10*binarray.shape[1]/100) == 0: print('doing time step ' + str(i) + 'of' + str(binarray.shape[1]) + 'for fish ') #+ str(y))
        i = i+1
        cid = np.where(binarray[:,t] > 0)[0]  #cid = cells active at current time point
    
            
        #LOOP THROUGH EACH ACTIVE CELL
        #-------------------------------
        #-------------------------------
        for c in cid:            #loop through all active cells at this time point

            if pkg[c,t] == 0:    #only find non-marked cells
                if len(np.intersect1d(np.where(nnbarray[c,:] > 0)[0], cid) > 2): #if >2 neighbours active
                    marker = marker + 1  
                    pkg[c,t] = marker  #mark active non-marked cell with new marker value
                       

            #LOCATE ALL NEIGHBOURS
            #----------------------------
            #----------------------------
            neighbour = np.where(nnbarray[c,:] > 0)[0]  #return indeces of current cell neighbours
            neighbouron  = np.intersect1d(cid,neighbour) #indeces of active cells in t, and also neighbours of c
            where0 = np.where(pkg[neighbouron,t] == 0)[0] #neighbours not already part of an avalanche
                
            #CONVERT NEIGHBOURS WHO ARE ALREADY PART OF AN AVALANCHE
            #-------------------------------------------------------
            #-------------------------------------------------------

            if len(where0) < len(neighbouron): #if any cells are already part of another avalanche
                oldav = np.unique(pkg[neighbouron, t]) #all avalanche values from neighbours
                firstav = np.min(oldav[np.where(oldav > 0)])   #minimum avalanche value that is not 0
                    
                    #define which cells we want to combine
                realav =  oldav[np.where(oldav > 0)] #all avalanche values that are not 0
                uniteav = np.where(pkg[:,t]==realav[:,None])[1] #indeces of all cells that need to be connected
                pkg[uniteav,t] = firstav #convert all current cell neighbours and their active neighbours 
                pkg[c,t] = firstav #also convert current cell
                    
                #GO BACK IN TIME AND CONVERT
                #----------------------------
                #----------------------------
                convertav = realav[1:] #avalanche numbers needing to be converted
                if t < 30:
                    time = t-1
                
                if t > 29:
                    time = 30
                        
                for e in range(convertav.shape[0]):
                    for timemachine in range(1, time): #loop through max possible time of previous avalanche
                        fill = np.where(pkg[:,t-timemachine] == convertav[e])[0]
                        if fill.shape[0] > 0:
                            pkg[fill, t-timemachine] = firstav 
                                    
            #CONVERT NEIGHBOURS WHO ARE NOT PART OF AN AVALANCHE
            #-------------------------------------------------------
            #-------------------------------------------------------
            if len(where0) == len(neighbouron): #if all cells are not part of an avalanche
                pkg[neighbouron[where0],t] = pkg[c,t]  

            
        #SEE IF AVALANCHE CAN PROPAGATE TO NEXT TIME FRAME
        #-------------------------------------------------------
        #-------------------------------------------------------
        n_av = np.unique(pkg[:,t])  #returns the marker values for each avalanche at this time point
    
        for n in n_av: #loop through each avalanche in this time point
            if n > 0:
                cgroup = np.where(pkg[:,t] == n)[0] #cells that are in same avalanche at t
                cid2 = np.where(binarray[:,t+1] > 0) #cells in next time point that are active
                intersect = np.intersect1d(cgroup, cid2) #check if any of the same cells are active in next time point
                wherealso0 = np.where(pkg[intersect,t+1] == 0)[0] #here we find all cells that are active in both time frames, and that are not already part of another avalanche - and mark them as current avalanche
                pkg[intersect[wherealso0], t+1] = pkg[cgroup[0],t] #carry over value to next frame for those cells
      
    allmark = np.unique(pkg)[1:] #all unique marker values

    #CALCULATE AVALANCHE SIZE
    #-------------------------------------------------------
    #-------------------------------------------------------
    avsize = np.unique(pkg, return_counts = True)[1][1:] #return counts for each unique avalanche
    frameslist = np.zeros(avsize.shape[0]) #create empty frames list of same length

    #CALCULATE AVALANCHE DURATION
    #-------------------------------------------------------
    #-------------------------------------------------------
    avpertimelist = list(range(pkg.shape[1])) #empty list of length time frames

    for e in range(pkg.shape[1]): #loop through each time point in pkg
            avpertime = np.unique(pkg[:,e]) #unique marker value in each time point
            avpertimelist[e] = avpertime #fill list of unique values in each time point
                          
    #link entire recording together
    #-----------------------------------------------------------
    linktime = list(itertools.chain(*avpertimelist)) #vector of all unique marker values in each time bin linked together
    framesvec = np.unique(linktime, return_counts = True)[1][1:] #vector of number of frames for each consecutive avalanche

    #COMBINE AV SIZE AND DURATION INTO ONE ARRAY
    #-------------------------------------------------------
    #-------------------------------------------------------
    avsizecut = avsize[avsize >= 3]  #only select avalanches above 2
    avframescut = framesvec[avsize >=3]
    av = np.vstack((avsizecut, avframescut))      
    return(av, pkg)


#=======================================================================
def avalanche(nnb, bind, savepath,experiment): # duration = yes convergence (no back propagation, earliest avalanche consumes meeting avalanche, and later avalanche terminates), cells in t must be active in t+1)
#=======================================================================
    import numpy as np
    import os
    import itertools

#Calculate avalanche size + duration
#-----------------------------------
    binarray, oldav, firstav, realav, timemachine, convertav, fill, time = [],[],[],[],[],[],[],[]
    
    #LOOP THROUGH EACH FISH
    #---------------------------------
    #---------------------------------
    binarray, nnbarray, pkg = np.load(bind),np.load(nnb), np.zeros(np.load(bind).shape)
    i, marker, avcount = 0,0,0
        
    #LOOP THROUGH EACH TIME POINT
    #------------------------------
    #------------------------------
    for t in range(binarray.shape[1]-1): #loop through all time points
        if i% round(10*binarray.shape[1]/100) == 0: print('doing time step ' + str(i) + 'of' + str(binarray.shape[1]) + 'for fish ') #+ str(y))
        i = i+1
        cid = np.where(binarray[:,t] > 0)[0]  #cid = cells active at current time point
    
            
        #LOOP THROUGH EACH ACTIVE CELL
        #-------------------------------
        #-------------------------------
        for c in cid:            #loop through all active cells at this time point

            if pkg[c,t] == 0:    #only find non-marked cells
                if len(np.intersect1d(np.where(nnbarray[c,:] > 0)[0], cid) > 2): #if >2 neighbours active
                    marker = marker + 1  
                    pkg[c,t] = marker  #mark active non-marked cell with new marker value
                       

            #LOCATE ALL NEIGHBOURS
            #----------------------------
            #----------------------------
            neighbour = np.where(nnbarray[c,:] > 0)[0]  #return indeces of current cell neighbours
            neighbouron  = np.intersect1d(cid,neighbour) #indeces of active cells in t, and also neighbours of c
            where0 = np.where(pkg[neighbouron,t] == 0)[0] #neighbours not already part of an avalanche
                
            #CONVERT NEIGHBOURS WHO ARE ALREADY PART OF AN AVALANCHE
            #-------------------------------------------------------
            #-------------------------------------------------------

            if len(where0) < len(neighbouron): #if any cells are already part of another avalanche
                oldav = np.unique(pkg[neighbouron, t]) #all avalanche values from neighbours
                firstav = np.min(oldav[np.where(oldav > 0)])   #minimum avalanche value that is not 0
                    
                #define which cells we want to combine
                realav =  oldav[np.where(oldav > 0)] #all avalanche values that are not 0
                uniteav = np.where(pkg[:,t]==realav[:,None])[1] #indeces of all cells that need to be connected
                pkg[uniteav,t] = firstav #convert all current cell neighbours and their active neighbours 
                pkg[c,t] = firstav #also convert current cell
                    
                #GO BACK IN TIME AND CONVERT
                #----------------------------
                #----------------------------
                convertav = realav[1:] #avalanche numbers needing to be converted
                if t < 30:
                    time = t-1
                
                if t > 29:
                    time = 30
                        
                for e in range(convertav.shape[0]):
                    for timemachine in range(1, time): #loop through max possible time of previous avalanche
                        fill = np.where(pkg[:,t-timemachine] == convertav[e])[0]
                        if fill.shape[0] > 0:
                            pkg[fill, t-timemachine] = firstav 
                                    
            #CONVERT NEIGHBOURS WHO ARE NOT PART OF AN AVALANCHE
            #-------------------------------------------------------
            #-------------------------------------------------------
            if len(where0) == len(neighbouron): #if all cells are not part of an avalanche
                pkg[neighbouron[where0],t] = pkg[c,t]  

            
        #SEE IF AVALANCHE CAN PROPAGATE TO NEXT TIME FRAME
        #-------------------------------------------------------
        #-------------------------------------------------------
        n_av = np.unique(pkg[:,t])  #returns the marker values for each avalanche at this time point
    
        for n in n_av: #loop through each avalanche in this time point
            if n > 0:
                cgroup = np.where(pkg[:,t] == n)[0] #cells that are in same avalanche at t
                cid2 = np.where(binarray[:,t+1] > 0) #cells in next time point that are active
                intersect = np.intersect1d(cgroup, cid2) #check if any of the same cells are active in next time point
                wherealso0 = np.where(pkg[intersect,t+1] == 0)[0] #here we find all cells that are active in both time frames, and that are not already part of another avalanche - and mark them as current avalanche
                pkg[intersect[wherealso0], t+1] = pkg[cgroup[0],t] #carry over value to next frame for those cells
      
    allmark = np.unique(pkg)[1:] #all unique marker values

    #CALCULATE AVALANCHE SIZE
    #-------------------------------------------------------
    #-------------------------------------------------------
    avsize = np.unique(pkg, return_counts = True)[1][1:] #return counts for each unique avalanche
    frameslist = np.zeros(avsize.shape[0]) #create empty frames list of same length

    #CALCULATE AVALANCHE DURATION
    #-------------------------------------------------------
    #-------------------------------------------------------
    avpertimelist = list(range(pkg.shape[1])) #empty list of length time frames

    for e in range(pkg.shape[1]): #loop through each time point in pkg
            avpertime = np.unique(pkg[:,e]) #unique marker value in each time point
            avpertimelist[e] = avpertime #fill list of unique values in each time point
                          
    #link entire recording together
    #-----------------------------------------------------------
    linktime = list(itertools.chain(*avpertimelist)) #vector of all unique marker values in each time bin linked together
    framesvec = np.unique(linktime, return_counts = True)[1][1:] #vector of number of frames for each consecutive avalanche

    #COMBINE AV SIZE AND DURATION INTO ONE ARRAY
    #-------------------------------------------------------
    #-------------------------------------------------------
    avsizecut = avsize[avsize >= 3]  #only select avalanches above 2
    avframescut = framesvec[[avsize >=3]]
    av = np.vstack((avsizecut, avframescut))
    return(av, pkg)


#==========================================================================
def powerfit(data, cutoff):
#==========================================================================
    import os
    import numpy as np
    import glob
    import powerlaw
    liklist = list(range(3)) 
    param = np.zeros((8))
    maxi = []
    #Truncated power law
    maxi = np.max(np.unique(data, return_counts = True)[0][np.unique(data, return_counts = True)[1] > cutoff]) #fit power law max - maximum value that appears more than 3 times
    fit = powerlaw.Fit(data, discrete = True, xmax = maxi) #initialise data for fit
    tp_alpha = fit.truncated_power_law.alpha
    tp_lam = fit.truncated_power_law.Lambda
    liklist[0] =  powerlaw.truncated_power_law_likelihoods(data, tp_alpha ,tp_lam, fit.xmin, xmax = maxi, discrete = True)
    param[0],param[1], param[2], param[3]= tp_alpha, tp_lam, fit.xmin, maxi
    if maxi - fit.xmin > 3:
        R, p = fit.distribution_compare('truncated_power_law', 'lognormal', normalized_ratio=True)
        param[4],param[5] = R,p
        R, p = fit.distribution_compare('truncated_power_law', 'exponential', normalized_ratio=True)
        param[6],param[7] = R,p

        #Lognormal
    fit = powerlaw.Fit(data, discrete = True) #initialise data for fit
    log_mu = fit.lognormal.mu
    log_sigma = fit.lognormal.sigma
    liklist[1] =  powerlaw.lognormal_likelihoods(data, log_mu ,log_sigma, fit.xmin, discrete = True)
    
        #Exponential
    fit = powerlaw.Fit(data, discrete = True) #initialise data for fit
    exp_lam = fit.exponential.Lambda
    liklist[2] =  powerlaw.exponential_likelihoods(data, exp_lam, fit.xmin, discrete = True)
    return(param, liklist)


#=======================================================================
def branch(pkg, av): # branching ratio calculation
#=======================================================================
    import numpy as np
    import os
    branchmean = []
    brancharr = np.zeros((np.int(np.max(pkg)), np.max(av[1])))
    i = 0
        
    for t in range(pkg.shape[1]): #loop through all time points
        if t == pkg.shape[1]-1:
            break
        n1 = np.unique(pkg[:,t])  #unique marker values at each time point
        n2 = np.unique(pkg[:,t+1]) 
        nx = np.intersect1d(n1, n2) #marker values that continue to next time frame
    
        if i% round(10*pkg.shape[1]/100) == 0: print('doing time step ' + str(i) + ' of ' + str(pkg.shape[1]))
        i = i+1

        for mark in nx[1:]: #loop through each marker value at this time point (only if marker active in next time point)
            if mark == brancharr.shape[0]:
                continue
            mark = np.int(mark)
            ancestor = np.unique(pkg[:,t], return_counts = True)[1][np.where(np.unique(pkg[:,t], return_counts = True)[0] == mark)[0]][0] #number of cells in that avalanche for that marker value at time point t  
            descend = np.unique(pkg[:,t+1], return_counts = True)[1][np.where(np.unique(pkg[:,t+1], return_counts = True)[0] == mark)[0]][0] #same as above for next time point
            brancharr[mark, np.where(brancharr[mark] == 0)[0][0]] = (descend/ancestor)
    branchmean = np.mean(brancharr[np.where(brancharr > 0)])
    return(branchmean)





#====================================================================================================
def ks_compare(distlist, mean_dist, bln_dist, choose, shape): #Ks distance groupwise statistical test
#=====================================================================================================
    import numpy as np
    from scipy import stats
    ks_p = np.zeros(2)
    ks = np.zeros(shape)
    for i in range(len(distlist)):
        ks[i] = stats.ks_2samp(mean_dist, np.load(distlist[i])[choose])[0]
    mean = np.mean(ks)
    sd = 1.96*np.std(ks)
    ks_diff = stats.ks_2samp(bln_dist[choose], mean_dist)[0]
    if ks_diff > mean + sd or ks_diff < mean - sd:
        p = 0
    else:
        p = 1
    ks_p = [ks_diff, p]
    return(ks_p)

#=====================================================================================================    
def euclidean_mat(input_coord, spatial_conversion):
#=====================================================================================================
    trans_coord = np.multiply(input_coord, spatial_conversion)
    mat = np.zeros((trans_coord.shape[0], trans_coord.shape[0]))
    for f in range(trans_coord.shape[0]):
        for x in range(trans_coord.shape[0]):
            mat[f,x] =  np.linalg.norm(trans_coord[f] - trans_coord[x])
            
    return(mat) 
        
#PLOT
#------------
#------------
        
#=======================================================================
def cellplot(Ftm, Fdrop, experiment, fnum, prefix, condition, plane, cell, xshift, yshift): # Plot cells and neighbours over image 
#=======================================================================
    import os
    import glob
    import numpy as np
    from matplotlib import pyplot as plt
    
    
    Freg = Ftm + 'Project/' + experiment + '-' + fnum + prefix
    os.chdir(Freg)
    opslist = sorted(glob.glob('*' + condition + '*' +'plane' + str(plane) + '*ops.npy'))
    os.chdir(Fdrop + 'Project/' + experiment)
    coordlist = sorted(glob.glob('*' + fnum +  '*' + condition + '*realcoord.npy'))
    nnblist = sorted(glob.glob('*' + fnum +  '*' + condition + '*nnb.npy'))
    
    if len(opslist) + len(coordlist) + len(nnblist) >3:
        print('More than one fish image loaded')
        
    # Plot
    #------------------------------------------------------
    ops = np.load(Freg + opslist[0])
    ops = ops[()]
    raw = ops['meanImg']

    # Pull out data from fish structure
    #----------------------------------------------------------------------------------------
    cs = np.load(coordlist[0])             # 3D array of xyz coordinates
    ci   = np.where(cs[:,2] == plane)[0]    # Index of plane coordinates in long list
    nnb = np.load(nnblist[0])

    # Actual plotting routines
    #----------------------------------------------------------------------------------------
    plt.figure(figsize = (15,15))
    plt.imshow(raw)
    plt.scatter(cs[ci,0]+xshift, cs[ci,1]-yshift, s = 10, c = nnb[ci[cell],ci], cmap = 'prism')
    plt.show()

            

            
#===============================================
def power_plot(data, data_inst, fig, units):
#===============================================
    from powerlaw import plot_pdf, Fit, pdf
    annotate_coord = (-.4, .95)
    ax1 = fig.add_subplot(n_graphs,n_data,data_inst)
    x, y = pdf(data, linear_bins=True)
    ind = y>0
    y = y[ind]
    x = x[:-1]
    x = x[ind]
    ax1.scatter(x, y, color='r', s=.5)
    plot_pdf(data[data>0], ax=ax1, color='b', linewidth=2)
    from pylab import setp
    setp( ax1.get_xticklabels(), visible=False)

    if data_inst==1:
        ax1.annotate("A", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)

    
    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    ax1in = inset_axes(ax1, width = "30%", height = "30%", loc=3)
    ax1in.hist(data, normed=True, color='b')
    ax1in.set_xticks([])
    ax1in.set_yticks([])

    
    ax2 = fig.add_subplot(n_graphs,n_data,n_data+data_inst, sharex=ax1)
    plot_pdf(data, ax=ax2, color='b', linewidth=2)
    fit = Fit(data, xmin=1, discrete=True)
    fit.power_law.plot_pdf(ax=ax2, linestyle=':', color='g')
    p = fit.power_law.pdf()

    ax2.set_xlim(ax1.get_xlim())
    
    fit = Fit(data, discrete=True)
    fit.power_law.plot_pdf(ax=ax2, linestyle='--', color='g')
    from pylab import setp
    setp( ax2.get_xticklabels(), visible=False)

    if data_inst==1:
       ax2.annotate("B", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)        
       ax2.set_ylabel(u"p(X)")# (10^n)")
        
    ax3 = fig.add_subplot(n_graphs,n_data,n_data*2+data_inst)#, sharex=ax1)#, sharey=ax2)
    fit.power_law.plot_pdf(ax=ax3, linestyle='--', color='g')
    fit.exponential.plot_pdf(ax=ax3, linestyle='--', color='r')
    fit.plot_pdf(ax=ax3, color='b', linewidth=2)
    
    ax3.set_ylim(ax2.get_ylim())
    ax3.set_xlim(ax1.get_xlim())
    
    if data_inst==1:
        ax3.annotate("C", annotate_coord, xycoords="axes fraction", fontproperties=panel_label_font)

    ax3.set_xlabel(units)
    
    
def marglik_power_loglik(data, npart):
    import numpy as np
    #Size
    sizes=data[0,:]
    M=len(sizes)
    a=min(sizes) #define xmin
    b=max(sizes) #define xmax
    size_ln=isfn.IS_LN(npart, sizes, M, a, b)
    size_po=isfn.IS(npart, sizes, M, a, b)
    
    #Dur
    sizes=data[1,:]
    a=2 #define xmin
    b=max(sizes) #define xmax
    M=len(sizes[np.where(sizes>a-1)])
    dur_ln=isfn.IS_LN(npart, sizes, M, a, b)
    dur_po=isfn.IS(npart, sizes, M, a, b)
    return(size_po, dur_po, size_ln, dur_ln)

def DCC(av):
    from matplotlib import pyplot as plt
    import numpy as np
    av_size = av[0]
    av_dur = av[1]
    ml = marglik_power_loglik(av, 2000)
    size_e = ml[0][0]
    dur_e = ml[1][0]
    fig, axarr = plt.subplots(figsize = (7,5))
    av_size = av_size
    av_dur = (1/2.73)*av_dur

    size_vec, dur_vec = [],[]
    for e in np.unique(av_dur):
        size_vec = np.append(size_vec, np.mean(av_size[np.where(av_dur == e)])) 
        dur_vec = np.append(dur_vec, e)

    xaxis = np.unique(dur_vec)
    yaxis = size_vec
    axarr.plot(xaxis[:len(xaxis)-1], yaxis[:len(xaxis)-1], '-', linewidth = 1.5, alpha = 1)
    fit_beta,c = np.polyfit(np.log10(xaxis[:len(xaxis)-1]), np.log10(yaxis[:len(xaxis)-1]), 1)
    plt.close(fig)
    
    pred_beta = (dur_e - 1)/(size_e - 1)
    dcc = np.abs(fit_beta - pred_beta)
    return(dcc)