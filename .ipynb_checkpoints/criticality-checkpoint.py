#PROCESS
#------------
#------------

#=======================================================================
def order(datalist, n, cond, mode, block): # Select which fish data to visualise
#=======================================================================
    import numpy as np
    
    # Plot longest block of data
    #-------------------------------------------------------------
    if mode == 'longest':
        mylist = datalist
        sublist = list(range(n*cond))
        count = 0

        for i in range(len(mylist)):
            if i == len(mylist)-1:
                sublist[count] = mylist[i]
            else:
                name1 = mylist[i][:mylist[i].find('run')-1]
                name2 = mylist[i+1][:mylist[i+1].find('run') -1]
                if name1 != name2:
                    sublist[count] = mylist[i]
                    count+=1
        return(sublist)
    
        # Plot all 30 minute block of data
        #-------------------------------------------------------------
    if mode == 'half':
        mylist =datalist
        sublist = list(range((n*cond)))
        count = 0
        for i in range(len(mylist)):
            if '30.npy' in mylist[i]: 
                sublist[count] = mylist[i]
                count+=1
            if '30' + str(block) + '.npy' in mylist[i]: 
                sublist[count] = mylist[i]
                count+=1
        return(sublist)
    

#=======================================================================
def neighbour(Fdrop, experiment, coordlist, rng, dim, cnt): # Select which fish data to visualise
#=======================================================================
    import numpy as np
    import os
    nnblist = list(range(len(coordlist)))   #nearest neighbour list containing nnb for each fish
    #Loop through all fish
    #----------------------
    for i in range(len(coordlist)):
        coord = np.load(coordlist[i])[:,:3]
        
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
                print("Doing row " + str(r) + " of " + str(coord.shape[0]) + " for " + coordlist[i][:coordlist[i].find('sess')-9] + '_' + coordlist[i][coordlist[i].find('dpf')+4:coordlist[i].find('run')-1])
            
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
        nnblist[i] = nnb
        np.save(Fdrop + 'Project/' + experiment +os.sep + coordlist[i][:coordlist[i].find('run')+6] + '_' + 'nnb.npy', nnb)
    return(nnblist)


#=======================================================================
def hist(Fdrop, experiment, sublist, mode): # Select which fish data to visualise
#=======================================================================
    import numpy as np
    import os
    
    if mode == 'size':
        blnlist = list(range(np.int(len(sublist)/3)))
        counter = 0
        for i in range(len(sublist)):
            if 'BLN' in sublist[i]:
                blnlist[counter] = sublist[i]
                counter +=1

        # Define histogram parameters
        #-----------------------------------------------------------
        countlist = list(range(len(sublist)))
        countblnlist = list(range(len(sublist)))
        htlist = list(range(len(sublist)))
        maxlist = list(range(len(sublist)))
        minlist = list(range(len(sublist)))

        #binning for ptz conditions
        for i in range(len(sublist)):
            countlist[i] = np.load(sublist[i])
            maxlist = max(map(lambda x: x, countlist[i]))
            minlist = min(map(lambda x: x, countlist[i]))
        maxi = np.max(maxlist)
        mini = np.min(minlist)
        bind = np.linspace(mini, maxi, 500)

        #binning for bln conditions
        for i in range(len(blnlist)):
            countblnlist[i] = np.load(blnlist[i])
            blnmaxlist = max(map(lambda x: x, countblnlist[i]))
            blnminlist = min(map(lambda x: x, countblnlist[i]))
        blnmaxi = np.max(blnmaxlist)
        blnmini = np.min(blnminlist)
        bindbln = np.linspace(blnmini, blnmaxi, 500)

        # Make histogram bins
        #-----------------------------------------------------------
        for y in range(len(htlist)):
            avdist = countlist[y]
            avdistcut  = avdist[avdist >= 2]
            if 'BLN' in sublist[y]:
                hist = np.histogram(avdistcut, bins = bindbln)
                htlist[y] = hist
                np.save(Fdrop + 'Project/' + experiment +os.sep + sublist[y][:sublist[y].find('run')+6] + '_' + 'avsizehist.npy', hist)
            else:
                hist = np.histogram(avdistcut, bins = bind)
                htlist[y] = hist
                np.save(Fdrop + 'Project/' + experiment +os.sep + sublist[y][:sublist[y].find('run')+6] + '_' + 'avsizehist.npy', hist)
        return(htlist)
        
    if mode == 'dur':     
        countlist = list(range(len(sublist)))
        htlist = list(range(len(sublist)))
        durmaxlist = list(range(len(sublist)))

        for i in range(len(sublist)):
            countlist[i] = np.load(sublist[i])
            durmaxlist[i] = len(countlist[i][np.where(countlist[i] > 0) [0]])
            durmaxi = np.max(durmaxlist)
            xlist = np.linspace(0, durmaxi, durmaxi+1)

        for i in range(len(sublist)):
            htlist[i] = countlist[i][:len(xlist)]
            np.save(Fdrop + 'Project/' + experiment +os.sep + sublist[i][:sublist[i].find('run')+6] + '_' + 'avdurhist.npy', htlist[i])
        return(htlist)

#ANALYSIS
#------------
#------------
#=======================================================================
def avsize(nnblist, binlist, Fdrop, experiment): # calculate avalanche sizes, size = each distinct avalanche in a distinct time frame is counted
#=======================================================================
    import numpy as np
    import os 

#Calculate avalanche size + duration
#-----------------------------------
    avlist = list(range(len(nnblist)))   
    pkglist = list(range(len(nnblist))) 
    
    #Loop through all fish
    #-----------------------------------------------------------------
    for y in range(len(avlist)):
        binarray = np.load(binlist[y])
        nnbarray = np.load(nnblist[y])
        pkg    = np.zeros(binarray.shape) #peak groups by timebin
        act_t  = np.where(np.sum(binarray, axis=0) > 3)[0] #Time points with at least 3 cells active
        i = 0 
        for t in range(binarray.shape[1]-1): #loop through all time points
            if i% round(10*binarray.shape[1]/100) == 0: print('doing time step ' + str(i) + ' of ' + str(binarray.shape[1]) + 'for fish ' + str(y))
            i = i+1    
            
            #Label each time point (of >3 cells active) with marker value and add to it as it grows
            #--------------------------------------------------------------------------------------
            marker  = 1    #set a marker
            cid = np.where(binarray[:,t] > 0)[0] #cid is list of cells at this time point (>3 cells active) that are firing at time frame t
            
            #Loop through all cells in current time frame with a cell firing
            #if currently unlabelled, label with marker + 1
            #fill empty matrix with ones as all starting points
            #--------------------------------------------------------------------------------------        
            for c in cid:  
                if pkg[c,t] == 0:       
                    marker = marker + 1
                    pkg[c,t] = marker    
            
            #Find all neighbours
            #nb = indeces of neighbours of cell c
            #tgr = indeces that are common to both cid and nb - cells that are active and neighbours 
            #-------------------------------------------------------------------
                nb   = np.where(nnbarray[c,:] > 0)[0]  
                tgr  = np.intersect1d(cid,nb)   
                
                #Fill all pkg at time point t with marker value so they are the same
                #-------------------------------------------------------------------
                pkg[tgr,t] = pkg[c,t]   
        # For each time point count number of coactivations
        #----------------------------------------------------------------------------
        print('Now calculating avalanche size')
        if binarray.shape[1] == 4914:
            av = np.array([])
            for t in range(binarray.shape[1]): 
                comms = np.unique(pkg[:,t][pkg[:,t] > 0])   
                for c in comms:
                    av = np.append(av, pkg[:,t][pkg[:,t] == c].shape[0])
            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')] + 'avsizelist30.npy', av)
            
        if binarray.shape[1] == 9828:
            av = np.array([])
            av1 = np.array([])
            av2 = np.array([])
            #loop through all time points
            #find all unique indeces at each time point where there is a value of 1 for a peak
            #append these coactivation values together to calculate total number of activations per time point
            #--------------------------------------------------------------------------------
            for t in range(binarray.shape[1]): 
                comms = np.unique(pkg[:,t][pkg[:,t] > 0])   
                for c in comms:
                    av = np.append(av, pkg[:,t][pkg[:,t] == c].shape[0]) 
                    if t < np.int(binarray.shape[1]/2):
                        av1 = np.append(av1, pkg[:,t][pkg[:,t] == c].shape[0])
                    if t > (np.int(binarray.shape[1]/2) - 1):
                        av2 = np.append(av2, pkg[:,t][pkg[:,t] == c].shape[0])
            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')] + 'avsizelist60.npy', av)
            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')] + 'avsizelist30a.npy', av1)
            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')] + 'avsizelist30b.npy', av2)
        pkglist[y] = pkg
        avlist[y] = av
        np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')] + 'avsizepkg.npy', pkg)
    return(avlist, pkglist)

#=======================================================================
def avduration(nnblist, binlist, Fdrop, experiment): # calculate avalanche duration, duration = normal (no convergence, cells in t must be active in t+1)
#=======================================================================
    import numpy as np
    import os
    import itertools

#Calculate avalanche size + duration
#-----------------------------------
    avlist = list(range(len(nnblist)))   
    pkglist = list(range(len(nnblist))) 
    for y in range(len(avlist)):
        binarray = np.load(binlist[y])
        nnbarray = np.load(nnblist[y])
        pkg    = np.zeros(binarray.shape) #peak groups by timebin
        i = 0
        marker = 0
        
        for t in range(binarray.shape[1]-1): #loop through all time points
            if i% round(10*binarray.shape[1]/100) == 0: print('doing time step ' + str(i) + ' of ' + str(binarray.shape[1]) + ' for fish ' + str(y))
            i = i+1
            cid = np.where(binarray[:,t] > 0)[0]  #cid = cells active at current time point
    
    #mark all avalanches in time point with different number
    #-------------------------------------------------------------------------------------
            for c in cid:            #loop through all active cells at this time point
                                     #mark all cells and its neighbours with same value 
                if pkg[c,t] == 0:
                    if len(np.intersect1d(np.where(nnbarray[c,:] > 0)[0], cid) > 2): #if >2 neighbours active
                        marker = marker + 1
                        pkg[c,t] = marker
        
        # Find all neighbours
        #-------------------------------------------------------------------
                neighbour = np.where(nnbarray[c,:] > 0)[0]  #return indeces of current cells neighbours
                neighbouron  = np.intersect1d(cid,neighbour)    #indeces of active cells in t, and also neighbours of c
                where0 = np.where(pkg[neighbouron,t] == 0)[0]
                pkg[neighbouron[where0],t] = pkg[c,t]  
    
    #mark all continuing avalanche with marker number in next frame
    #-------------------------------------------------------------------------------------
            n_av = np.unique(pkg[:,t])  #returns the marker values for each avalanche at this time point
            for n in n_av: #loop through each avalanche in this time point
                if n > 0:
                    cgroup = np.where(pkg[:,t] == n)[0] #cells that are in same avalanche at t
                    cid2 = np.where(binarray[:,t+1] > 0) #cells in next time point that are active
                    intersect = np.intersect1d(cgroup, cid2) #check if any of the same cells are active in next time point
                    wherealso0 = np.where(pkg[intersect,t+1] == 0)[0] #here we find all cells that are active in both time frames, and that are not part of an avalanche - and mark them as avalanche
                    pkg[intersect[wherealso0], t+1] = pkg[cgroup[0],t] #carry over value to next frame for those cells
        
        #Calculate unique avalanche marker values for each time point
        #------------------------------------------------------------
        print('Now calculating avalanche duration')
        if binarray.shape[1] == 4914:
            uniqvalist = list(range(pkg.shape[1])) #empty list of length time frames
            for e in range(pkg.shape[1]): #loop through each time point in pkg
                uniqval = np.unique(pkg[:,e]) #unique marker value in each time point
                uniqvalist[e] = uniqval #fill list of unique values in each time point
            
            #link entire recording together
            #-----------------------------------------------------------
            linktime = list(itertools.chain(*uniqvalist)) #vector of all unique marker values in each time bin linked together
            
            #find the max length of avalanches
            #----------------------------------
            noavlist = 0
            index = []
            for f in range(pkg.shape[1]): #loop through all time points in pkg to find longest length
                if f > 0:  
                    if len(np.where(np.unique(linktime, return_counts = 'True')[1] == f)[0]) == 0: #find how many avalanches occur for each number of frames - each time you find an avalanche duration that has 0 avalanches - add to z list - keep iterating however as some avalanche lengths may skip frames - e.g. you may have no avalanche 5 frames long, but some 6.. etc long
                        noavlist +=1
                        if noavlist > 50: #once you have reached 50 empty avalanche durations set vector length to this and break
                            index = f 
                            break
        
        #calculate how long each marker value is repeated for consecutively - any repeat must be a consecutive time frame
            time = np.zeros(index, dtype = 'int') #vector of length index
            for o in range(index):
                time[o] = len(np.where(np.unique(linktime, return_counts = 'True')[1] == o)[0])
            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')] + 'avdurlist30.npy', time)
    
    
        if binarray.shape[1] == 9828:
               
            uniqvalist = list(range(pkg.shape[1])) #empty list of length time frames
            uniqvalist1 = list(range(np.int(pkg.shape[1]/2))) #empty list of length time frames
            uniqvalist2 = list(range(np.int(pkg.shape[1]/2))) #empty list of length time frames
    
            for e in range(pkg.shape[1]): #loop through each time point in pkg
                uniqval = np.unique(pkg[:,e]) #unique marker value in each time point
                uniqvalist[e] = uniqval #fill list of unique values in each time point
                uniqvalist1 = uniqvalist[:np.int(pkg.shape[1]/2)]
                uniqvalist2 = uniqvalist[np.int(pkg.shape[1]/2):]
                
                #link entire recording together
                #-----------------------------------------------------------
            linktime = list(itertools.chain(*uniqvalist)) #vector of all unique marker values in each time bin linked together
            linktime1 = list(itertools.chain(*uniqvalist1)) 
            linktime2 = list(itertools.chain(*uniqvalist2)) 
                #find the max length of avalanches
                #----------------------------------
            noavlist = 0
            index = []
            for f in range(pkg.shape[1]): #loop through all time points in pkg to find longest length
                if f > 0:  
                    if len(np.where(np.unique(linktime, return_counts = 'True')[1] == f)[0]) == 0: #find how many avalanches occur for each number of frames - each time you find an avalanche duration that has 0 avalanches - add to z list - keep iterating however as some avalanche lengths may skip frames - e.g. you may have no avalanche 5 frames long, but some 6.. etc long
                        noavlist +=1
                        if noavlist > 50: #once you have reached 50 empty avalanche durations set vector length to this and break
                            index = f 
                            break
        #calculate how long each marker value is repeated for consecutively - any repeat must be a consecutive time frame
            time = np.zeros(index, dtype = 'int')
            time1 = np.zeros(index, dtype = 'int') 
            time2 = np.zeros(index, dtype = 'int') #vector of length index
            for o in range(index):
                time[o] = len(np.where(np.unique(linktime, return_counts = 'True')[1] == o)[0])
                time1[o] = len(np.where(np.unique(linktime1, return_counts = 'True')[1] == o)[0])
                time2[o] = len(np.where(np.unique(linktime2, return_counts = 'True')[1] == o)[0])

            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')]  + 'avdurlist60.npy', time)
            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')]  + 'avdurlist30a.npy', time1)
            np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')]  + 'avdurlist30b.npy', time2)
        
                    
        avlist[y] = time
        pkglist[y] = pkg
        np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[y][:nnblist[y].find('nnb')] + 'avdurpkg.npy', pkg)
    return(avlist, pkglist)


        
#=======================================================================
def branch(pkglist, subdurlist, Fdrop, experiment): # calculate branching ratio
#=======================================================================
    import numpy as np
    import os
    
    branchlist = list(range(len(pkglist)))   
    for y in range(len(pkglist)):
        pkg = np.load(Fdrop + 'Project/' + experiment +  '/criticality/pkg_dur/' + os.sep + pkglist[y])
        dur = np.load(Fdrop + 'Project/' + experiment + '/criticality/duration/' +  os.sep + subdurlist[y])
    
        brancharr = np.zeros((np.int(np.max(pkg)), (np.where(dur == 0)[0][1])+100)) #empty array of size: number of avalanches x max number of frames
        list(range(np.int(np.max(pkg)))) #empty list - of dimensions - number of marker values x max duration of avalanche    
        i = 0
        for t in range(pkg.shape[1]): #loop through all time points
            if t == pkg.shape[1]-1:
                break
            n1 = np.unique(pkg[:,t])  #unique marker values at each time point
            n2 = np.unique(pkg[:,t+1])
            nx = np.intersect1d(n1, n2) #marker values that continue to next time frame
    
            if i% round(10*pkg.shape[1]/100) == 0: print('doing time step ' + str(i) + ' of ' + str(pkg.shape[1]) + ' for fish ' + str(y))
            i = i+1

            for mark in nx[1:]: #loop through each marker value at this time point (only if marker active in next time point)
                mark = np.int(mark)
                ancestor = np.unique(pkg[:,t], return_counts = True)[1][np.where(np.unique(pkg[:,t], return_counts = True)[0] == mark)[0]][0] #number of cells in that avalanche for that marker value at time point t  
                descend = np.unique(pkg[:,t+1], return_counts = True)[1][np.where(np.unique(pkg[:,t+1], return_counts = True)[0] == mark)[0]][0] #same as above for next time point
                brancharr[mark, np.where(brancharr[mark] == 0)[0][0]] = (descend/ancestor)
        branchlist[y] = np.mean(brancharr[np.where(brancharr > 0)])
        np.save(Fdrop + 'Project/' + experiment + os.sep + pkglist[y][:pkglist[y].find('run')+7] + 'branch.npy', branchlist[y])
    return(branchlist)


#=======================================================================
def exp(Fdrop, experiment, histlist, mode): # calculate critical exponent
#=======================================================================
    import numpy as np
    import os
    from scipy import stats

    ylist = list(range(len(histlist)))
    xlist =list(range(len(histlist)))
    linelist = list(range(len(histlist)))
    slopelist = list(range(len(histlist)))
    
    for i in range(len(histlist)):
        if mode == 'size':
            y0 = np.log(np.load(histlist[i])[0])
            x0 = np.log(np.load(histlist[i])[1][:-1])
        if mode == 'dur':
            y0 = np.log(np.load(histlist[i]))
            x0 = np.log(np.linspace(0, len(y0), len(y0))) 
        finiteymask = np.isfinite(y0)
        yclean = y0[finiteymask]
        xclean = x0[finiteymask]
        slope, intercept, r_value, p_value, std_err = stats.linregress(xclean,yclean)
        line0 = slope*x0+intercept
        ylist[i] = y0
        xlist[i] = x0
        linelist[i] = line0
        slopelist[i] = slope
        
        np.save(Fdrop + 'Project/' + experiment + os.sep + histlist[i][:histlist[i].find('run')+7] + mode + 'exponent.npy', slope)

        
#Goodness of fit
#Loglikelihood ratios - comparison of power law fit, to fit to other distributions
#Loglikelihood is faster, and more accurate than bootstrapping 
#(bootstrapping could find a power law with sufficient likelihood, but there may yet be better distribution untested)
#Hard to say if a distribution is really a power law (impossible to follow theoretical) 
#Can instead ask if it the best description available
#exponential distribution - minimum alternative, heavy tailed cannot be exponentially bounded
#fit object - list of supported distributions
#R = loglikelihood ratio between 2 candidate distributions - positive = first distribution, negative = second distribution
#p = significance value for ratio
#normalised ratio - normalises R by its SD - used to calculate p
#-------------------------------------------------------------------------------------------------------------------
#==========================================================================
def loglik(Fdrop, experiment, distlist, dist1, dist2, normratio = True):
#==========================================================================
    import powerlaw
    import numpy as np
    import os
    Rlist = list(range(len(distlist)))
    plist = list(range(len(distlist)))
    for i in range(len(distlist)):
        data = np.load(distlist[i])
        fit = powerlaw.Fit(data)
        R, p = fit.distribution_compare(dist1, dist2, normalized_ratio = normratio)
        Rlist[i] = R
        plist[i] = p
        np.save(Fdrop + 'Project/' + experiment + os.sep + distlist[i][:distlist[i].find('run')+6] + '_' + 'loglik.npy', R)
        np.save(Fdrop + 'Project/' + experiment + os.sep + distlist[i][:distlist[i].find('run')+6] + '_' + 'loglikp.npy', p)
    return(Rlist,plist)   
        
        
        
        
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

    
#=======================================================================
def avplot(histlist, mode, plot, n): # Plot av distribution
#=======================================================================
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt


    if mode == 'size':
        if plot == 'single':
            n = n*3
            fig, ax = plt.subplots(figsize = (12,12))
            plt.scatter(np.log(np.load(histlist[n])[1][:-1]), np.log(np.load(histlist[n])[0]/np.log(np.sum(np.load(histlist[n])[0]))), label = 'baseline', s = 100)
            plt.scatter(np.log(np.load(histlist[n+1])[1][:-1]), np.log(np.load(histlist[n+1])[0]/np.log(np.sum(np.load(histlist[n+1])[0]))), label = 'ptz5mm', s=100)
            plt.scatter(np.log(np.load(histlist[n+2])[1][:-1]), np.log(np.load(histlist[n+2])[0]/np.log(np.sum(np.load(histlist[n+2])[0]))), label = 'ptz20mm', s = 100, c = 'r')
            ax.legend()
            plt.show()
            
        # Multi plot size and calculate exponent
        #--------------------------------------------------------------------------------------
        if plot == 'all':
            f, axarr = plt.subplots(4,3,sharey=True, sharex=True, figsize = (20,20))
            f.subplots_adjust(hspace=0)
            axarr = axarr.flatten()
            
            icounter = 0
            plt.subplots_adjust(wspace=0, hspace=0)
            for i in range(int(len(histlist)/3)):
                axarr[i].scatter(np.log(np.load(histlist[icounter])[1][:-1]), np.log(np.load(histlist[icounter])[0]), label = 'baseline', s = 8)
                axarr[i].scatter(np.log(np.load(histlist[icounter+1])[1][:-1]), np.log(np.load(histlist[icounter+1])[0]), label = 'ptz5mm', s=8)
                axarr[i].scatter(np.log(np.load(histlist[icounter+2])[1][:-1]), np.log(np.load(histlist[icounter+2])[0]), label = 'ptz20mm', s = 8, c = 'r')
                icounter = icounter + 3

            axarr[7].set_xlabel('Avalanche size', fontsize = 40)
            axarr[3].set_ylabel('Probability', fontsize = 40)
            plt.show()
            
            
    # Define histogram parameters
    #-----------------------------------------------------------           
    if mode == 'dur':
            
        if plot == 'single':
            n = n*3
            fig, ax = plt.subplots(figsize = (12,12))
            xlist = np.linspace(0, len(np.load(histlist[n])), len(np.load(histlist[n])))
            plt.scatter(np.log(xlist), np.log(np.load(histlist[n]))/np.log(np.sum(np.load(histlist[n]))), label = 'baseline', s=50)
            plt.scatter(np.log(xlist), np.log(np.load(histlist[n+1]))/np.log(np.sum(np.load(histlist[n+1]))), label = 'ptz5mm', s=50)
            plt.scatter(np.log(xlist), np.log(np.load(histlist[n+2]))/np.log(np.sum(np.load(histlist[n+2]))), label = 'ptz20mm', s=50, c = 'r')
            plt.legend(loc='upper left')
            
            
        # Multiple plot duration
        #---------------------------------------------------------------------------------------
        if plot == 'all':

            f, axarr = plt.subplots(4,3,sharey=True, sharex=True, figsize = (20,20))
            f.subplots_adjust(hspace=0)
            axarr = axarr.flatten()
            icounter = 0

            plt.subplots_adjust(wspace=0, hspace=0)
            xlist = np.linspace(0, len(np.load(histlist[n])), len(np.load(histlist[n])))

            for i in range(int(len(histlist)/3)):
                axarr[i].scatter(np.log(xlist), np.log(np.load(histlist[icounter]))/np.log(np.sum(np.load(histlist[icounter]))), label = 'baseline', s = 15)
                axarr[i].scatter(np.log(xlist), np.log(np.load(histlist[icounter+1]))/np.log(np.sum(np.load(histlist[icounter+1]))), label = 'ptz5mm', s=15)
                axarr[i].scatter(np.log(xlist), np.log(np.load(histlist[icounter+2]))/np.log(np.sum(np.load(histlist[icounter+2]))), label = 'ptz20mm', s = 15, c = 'r')
                axarr[i].legend(loc='upper left')
                icounter = icounter + 3

            axarr[7].set_xlabel('Avalanche duration', fontsize = 40)
            axarr[3].set_ylabel('Probability', fontsize = 40)
            plt.show()

            
            

#==========================================================================
def powerfit(Fdrop, experiment, distlist, mode):
#==========================================================================

    import pylab
    pylab.rcParams['xtick.major.pad']='8'
    pylab.rcParams['ytick.major.pad']='8'
    from matplotlib import rc
    rc('font', family='sans-serif')
    rc('font', size=10.0)
    rc('text', usetex=False)
    from matplotlib.font_manager import FontProperties
    panel_label_font = FontProperties().copy()
    panel_label_font.set_weight("bold")
    panel_label_font.set_size(12.0)
    panel_label_font.set_family("sans-serif")

        #PLOT - log binning (blue) vs linear binning (red)
        #calculate density functions - pdf, cdf and ccdf 
        #plot commands - plot_pdf, plot_cdf, and plot_ccdf, takes matplotlib keywordargs
        #-------------------------------------------------------------------------------
    if mode == 'loglin':
        data = distlist
        for i in range(len(data)):
            fig, ax = plt.subplots(figsize = (5,5)) 
            figPDF = powerlaw.plot_pdf(np.load(data[i]), color='b', label = 'log') #default log binning - increases likelihood of observing a range of                                                                         #values in the tail 
#of distribution where some high values may occur, rather than binning it into single bin
            powerlaw.plot_pdf(np.load(data[i]), linear_bins=True, color='r', ax=figPDF)
####
            figPDF.set_ylabel("Probability", size = 20)
            figPDF.set_xlabel(r"Avalanche size", size = 20)
#PLOT - power law fit (dotted) to probability density (blue) and cumulative density (red)
#fit - powerlaw, lognormal
#-----------------------------------------------------------------------------------------

    if mode == 'fit':

        for i in range(len(praclist)):
            fig, ax = plt.subplots(figsize = (5,5))
            data = np.load(distlist[i])
            fit = powerlaw.Fit(data, discrete=True) #fit power law to your data
####
            figCCDF = fit.plot_pdf(color='b', linewidth=2) #pdf is better at visualising changes in the tail
            fit.power_law.plot_pdf(color='b', linestyle='--', ax=figCCDF) 
            fit.plot_ccdf(color='r', linewidth=2, ax=figCCDF) #ccdf does not require binning
            fit.power_law.plot_ccdf(color='r', linestyle='--', ax=figCCDF) 
####
            figCCDF.set_ylabel(u"p(X),  p(X≥x)")
            figCCDF.set_xlabel(r"Word Frequency")

            figname = 'FigCCDF'



            
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