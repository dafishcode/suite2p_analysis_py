#WHERE OLD CODE GOES TO DIE - AND BE STORED IN CASE
#=============================================


#============================================================================
#============================================================================
#AVALANCHES
#============================================================================
#============================================================================
#AVALANCHE DURATION - OLD WAY
#------------------------------

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
time = np.zeros(index, dtype = 'int')#vector of length index
for o in range(index):
    time[o] = len(np.where(np.unique(linktime, return_counts = 'True')[1] == o)[0])
        
                    
#avlist[y] = time
#pkglist[y] = pkg

#AVALANCHE SIZE - OLD WAY
#------------------------------

#Calculate size based on duration (each avalanche at each time point treated as a distinct event)
#-------------------------------------
newavlist = list(range(len(praclist)))   
for i in range(len(praclist)):
    pkg = np.load(pracpkg[i])

    binarray = np.load(binlist[i])
    av = np.array([])

            #loop through all time points
            #find all unique indeces at each time point where there is a value of 1 for a peak
            #append these coactivation values together to calculate total number of activations per time point
            #--------------------------------------------------------------------------------
    for t in range(binarray.shape[1]): 
        comms = np.unique(pkg[:,t][pkg[:,t] > 0])   
        for c in comms:
            av = np.append(av, pkg[:,t][pkg[:,t] == c].shape[0])    
    
    newavlist[i] = av
    np.save(Fdrop + 'Project/' + experiment + os.sep + nnblist[i][:nnblist[i].find('nnb')] + 'opracavsizelist.npy', newavlist[i])

#Old avalanche duration
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


#Old avalanche size calculation
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
#old histogram calculation
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
            countblnlist[i] = np.load(sublist[i])
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
                np.save(Fdrop + 'Project/' + experiment +os.sep + sublist[y][:sublist[y].find('run')+6] + '_' + 'opracavsizehist.npy', hist)
            else:
                hist = np.histogram(avdistcut, bins = bind)
                htlist[y] = hist
                np.save(Fdrop + 'Project/' + experiment +os.sep + sublist[y][:sublist[y].find('run')+6] + '_' + 'opracavsizehist.npy', hist)
        return(htlist)
        
    if mode == 'dur':     
        countlist = list(range(len(sublist)))
        htlist = list(range(len(sublist)))
        durmaxlist = list(range(len(sublist)))

        for i in range(len(sublist)):
            countlist[i] = np.load(sublist[i])
            durmaxlist[i] = np.max(np.where(countlist[i] > 0) [0])
        durmaxi = np.max(durmaxlist)
        xlist = np.linspace(0, durmaxi, durmaxi+3)

        for i in range(len(sublist)):
            addon = len(xlist) - countlist[i].shape[0]
            addzeros = np.zeros(addon)
            htlist[i] = np.append(countlist[i], addzeros)
            np.save(Fdrop + 'Project/' + experiment +os.sep + sublist[i][:sublist[i].find('run')+6] + '_' + 'opracavdurhist.npy', htlist[i])
        return(htlist)
    
    
# calculates number of cells in each distinct avalanche event (takes ages)
#=======================================================================
def distinctduration(nnblist, binlist): # duration = normal, size = calculate from unique cells in duration
#=======================================================================
    import numpy as np
    
#Calculate avalanche size + duration
#-----------------------------------
    avlist = list(range(len(nnblist)))
    distlist = list(range(len(nnblist)))
    
    for y in range(len(avlist)):
        binarray = np.load(binlist[y])
        nnbarray = np.load(nnblist[y])
        pkg    = np.zeros(binarray.shape) #peak groups by timebin
        i = 0
        marker = 0
        distinctlist = (np.zeros(1)) #list containing updated number of distinct cells per avalanche
        
        for t in range(binarray.shape[1]-1): #loop through all time points
            #if i% round(10*binarray.shape[1]/100) == 0: print('doing time step ' + str(i) + 'of' + str(binarray.shape[1]) + 'for fish ' + str(y))
            #i = i+1
            print('done time point' + ' ' + str(t))
            cid = np.where(binarray[:,t] > 0)[0]  #cid = cells active at current time point
    
    #mark all avalanches in time point with different number
    #-------------------------------------------------------------------------------------
            for c in cid:            #loop through all active cells at this time point
                                     #mark all cells and its neighbours with same value 
                if pkg[c,t] == 0:
                    if len(np.intersect1d(np.where(nnbarray[c,:] > 0)[0], cid) > 2): #if >2 neighbours active
                        marker = marker + 1
                        pkg[c,t] = marker
                        distinctlist = np.hstack((distinctlist, (np.zeros(1)))) #once each marker is created append empty element to distinctlist vector for future addition
                        
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
                    
                    #if avalanche ends at this time point then calculate unique cells in this avalanche
                    if np.sum(intersect) == 0:
                        distinctlist[(np.int(n))] = len(np.unique(np.where(pkg == n)[0]))                                     
        avlist[y] = pkg
        distlist[y] = distinctlist
        
    return(avlist, distinctlist)

#=======================================================================
def test2(nnblist, binlist): # duration = no convergence, cells in t need not be active in t+1, neighbours can be active
#=======================================================================
    import numpy as np
    
#Calculate avalanche size + duration
#-----------------------------------
    avlist = list(range(len(nnblist)))
    
    for y in range(len(avlist)):
        binarray = np.load(binlist[y])
        nnbarray = np.load(nnblist[y])
        pkg    = np.zeros(binarray.shape) #peak groups by timebin
        i = 0
        marker = 0
        
        for t in range(binarray.shape[1]-1): #loop through all time points
            if i% round(10*binarray.shape[1]/100) == 0: print('doing time step ' + str(i) + 'of' + str(binarray.shape[1]) + 'for fish ' + str(y))
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
                    wherealso0 = np.where(pkg[intersect,t+1] == 0)[0] #here we find all cells that are active in both time frames, and that are not already part of another avalanche - and mark them as current avalanche
                    
                    neighboursnow = np.where(nnbarray[cgroup,:] > 0)[1]   #all neighbours of currently active cells
                    intersectneighcellt1 = np.intersect1d(neighboursnow,cid2) #cells active in next time point and neighbours of current time frame
                    
                    pkg[intersect[wherealso0], t+1] = pkg[cgroup[0],t] #carry over value to next frame for cells in last time point who are also active in next time point
                    pkg[intersectneighcellt1,t+1] = pkg[cgroup[0],t] #carry over marker value for neighbours of cells active in last time point, who are active in next
          
        avlist[y] = pkg    
    return(avlist)

#kernel density estimate
#---------------
from scipy.stats import kde
trans = Fishcoordz[pracstack][0][:,:2].T
k = kde.gaussian_kde(trans)
xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))

#SPATIAL EXTENT
#------------------
#metrically scaled coordinates

dim = [.8, .8, 15]
coord = Fishcoordz6
mcs  = np.multiply(coord, dim)  


meandistance7 = np.zeros(130)
for i in range(130):
    #calculate mean coords for all avalanche
    distance = np.zeros(mcs[avcoord[i][0]].shape[0])
    meanx = np.mean(mcs[avcoord[i][0]][:,0])
    meany = np.mean(mcs[avcoord[i][0]][:,1])
    meanz = np.mean(mcs[avcoord[i][0]][:,2])
    meancoord = np.hstack((meanx,meany,meanz))
    #find euclidean distance for each coordinate
    #loop through each cell in avalanche|
    for j in range(mcs[avcoord[i][0]].shape[0]):
        distance[j] = np.linalg.norm(meancoord-mcs[avcoord[i][0]][j])
    meandistance7[i] = np.mean(distance)

    
    
#============================================================================
#============================================================================
#PLOT
#============================================================================
#============================================================================

#=======================================================================
def avplot(histlist, mode, plot, n): # Plot av distribution
#=======================================================================
    import numpy as np
    from scipy import stats
    import matplotlib.pyplot as plt
    import os


    if mode == 'size':
        if plot == 'single':
            n = n*3
            fig, ax = plt.subplots(figsize = (6,6))
            plt.scatter(np.log(np.load(histlist[n])[1][:-1]), np.log(np.load(histlist[n])[0]/np.sum(np.load(histlist[n])[0])), label = 'baseline', s = 20)
            plt.scatter(np.log(np.load(histlist[n+1])[1][:-1]), np.log(np.load(histlist[n+1])[0]/np.sum(np.load(histlist[n+1])[0])), label = 'ptz5mm', s=20)
            plt.scatter(np.log(np.load(histlist[n+2])[1][:-1]), np.log(np.load(histlist[n+2])[0]/np.sum(np.load(histlist[n+2])[0])), label = 'ptz20mm', s = 20, c = 'r')
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
                axarr[i].scatter(np.log(np.load(histlist[icounter])[1][:-1]), np.log(np.load(histlist[icounter])[0]/np.sum(np.load(histlist[icounter])[0])), label = 'baseline', s = 8)
                axarr[i].scatter(np.log(np.load(histlist[icounter+1])[1][:-1]), np.log(np.load(histlist[icounter+1])[0]/np.sum(np.load(histlist[icounter+1])[0])), label = 'ptz5mm', s=8)
                axarr[i].scatter(np.log(np.load(histlist[icounter+2])[1][:-1]), np.log(np.load(histlist[icounter+2])[0]/np.sum(np.load(histlist[icounter+2])[0])), label = 'ptz20mm', s = 8, c = 'r')
                icounter = icounter + 3

            axarr[7].set_xlabel('Avalanche size', fontsize = 40)
            axarr[3].set_ylabel('Probability', fontsize = 40)
            plt.show()
            
            
    # Define histogram parameters
    #-----------------------------------------------------------           
    if mode == 'dur':
            
        if plot == 'single':
            n = n*3
            fig, ax = plt.subplots(figsize = (6,6))
            xlist = np.linspace(0, len(np.load(histlist[n])), len(np.load(histlist[n])))
            plt.scatter(np.log(xlist), np.log(np.load(histlist[n]))/np.log(np.sum(np.load(histlist[n]))), label = 'baseline', s=50)
            plt.scatter(np.log(xlist), np.log(np.load(histlist[n+1]))/np.log(np.sum(np.load(histlist[n+1]))), label = 'ptz5mm', s=50)
            plt.scatter(np.log(xlist), np.log(np.load(histlist[n+2]))/np.log(np.sum(np.load(histlist[n+2]))), label = 'ptz20mm', s=50, c = 'r')
            plt.show()
            

            
            
            
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


#avalanche distribution plotting regimes
#-----------------------------------------


f, axarr = plt.subplots(3,3,sharey=True, sharex=True, figsize = (20,20))
f.subplots_adjust(hspace=0)
axarr = axarr.flatten()
icounter = 0

plt.subplots_adjust(wspace=0, hspace=0)

for i in range(int(len(critfold)/3)):
    axarr[i].scatter(np.log(ahlist[icounter][1][:-1]), np.log(ahlist[icounter][0])/np.log(np.sum(ahlist[icounter][0])), label = 'baseline', s = 50)
    axarr[i].scatter(np.log(ahlist[icounter+1][1][:-1]), np.log(ahlist[icounter+1][0])/np.log(np.sum(ahlist[icounter+1][0])), label = 'ptz5mm', s=50)
    axarr[i].scatter(np.log(ahlist[icounter+2][1][:-1]), np.log(ahlist[icounter+2][0])/np.log(np.sum(ahlist[icounter+2][0])), label = 'ptz20mm', s = 50, c = 'r')
    #axarr[i].legend(loc='upper left')
    icounter = icounter + 3



axarr[7].set_xlabel('Avalanche size', fontsize = 40)
axarr[3].set_ylabel('Probability', fontsize = 40)

#$10^1$','$10^2$', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$', '$10^8$', '$10^9$' ), size = 15, color = 'white')
#plt.yticks(np.arange(15), (

#axarr[8].set_xticklabels([(['1','$\mathregular{10^{-1}}$', '$\mathregular{10^{-0.7}}$', '$\mathregular{10^{-0.4}}$', '$\mathregular{10^{-0.2}}$','$\mathregular{10^{-0.09}}$', '$\mathregular{10^{0}}$'], fontsize = 15) 
#axarr[7].set_xticklabels(['-0.2', '-0.2','0', '0.2', '0.4', '0.6'], fontsize = 15)
axarr[6].set_xticklabels(['$\mathregular{2x10^2}}$', '$\mathregular{10^1}}$','$\mathregular{10^2}}$', '$\mathregular{10^3}}$', '$\mathregular{10^4}}$', '$\mathregular{2x10^2}}$'], fontsize = 15)
axarr[0].set_yticklabels(['1','$\mathregular{10^{-1}}$', '$\mathregular{10^{-0.7}}$', '$\mathregular{10^{-0.4}}$', '$\mathregular{10^{-0.2}}$','$\mathregular{10^{-0.09}}$', '$\mathregular{10^{0}}$'], fontsize = 15)
axarr[3].set_yticklabels(['1','$\mathregular{10^{-1}}$', '$\mathregular{10^{-0.7}}$', '$\mathregular{10^{-0.4}}$', '$\mathregular{10^{-0.2}}$','$\mathregular{10^{-0.09}}$', '$\mathregular{10^{0}}$'], fontsize = 15)
axarr[6].set_yticklabels(['1','$\mathregular{10^{-1}}$', '$\mathregular{10^{-0.7}}$', '$\mathregular{10^{-0.4}}$', '$\mathregular{10^{-0.2}}$','$\mathregular{10^{-0.09}}$', '$\mathregular{10^{0}}$'], fontsize = 15)
axarr[0].locator_params(axis='x', nbins=5)

axarr[0].yaxis.set_tick_params(labelsize = 20)
axarr[8].xaxis.set_tick_params(labelsize=20)
axarr[3].yaxis.set_tick_params(labelsize = 20)
axarr[6].xaxis.set_tick_params(labelsize=20)
axarr[6].yaxis.set_tick_params(labelsize = 20)
axarr[7].xaxis.set_tick_params(labelsize=20)
axarr[2].legend(loc = 1, markerscale = 2,prop={'size': 30})


axarr[8].plot(x0, line0, c = 'k', linestyle = '--', linewidth = 3)
axarr[8].plot(x1, line1, c = 'k', linestyle = '--', linewidth = 3)
axarr[8].plot(x2, line2, c = 'k', linestyle = '--', linewidth = 3)




f, axarr = plt.subplots(3,3,sharey=True, sharex=True, figsize = (20,20))
f.subplots_adjust(hspace=0)
axarr = axarr.flatten()

plt.subplots_adjust(wspace=0, hspace=0)

for i in range(int(len(datalist))):
    bln  = datalist[i].iloc[0,:]
    p5 = datalist[i].iloc[1,:]
    p20 = datalist[i].iloc[2,:]
    x = np.linspace(1, bln.shape[0], num = bln.shape[0])
    axarr[i].scatter(np.log(x), np.log(bln)/np.log(np.sum(bln)), label = 'baseline', s=50)
    axarr[i].scatter(np.log(x), np.log(p5)/np.log(np.sum(p5)), label = 'ptz5mm', s=50)
    axarr[i].scatter(np.log(x), np.log(p20)/np.log(np.sum(p20)), label = 'ptz20mm', s=50, c = 'r')
    #axarr[i].legend(loc='upper left')

axarr[7].set_xlabel('Avalanche duration', fontsize = 40)
axarr[3].set_ylabel('Probability', fontsize = 40)

axarr[6].set_xticklabels(['$\mathregular{2x10^2}}$', '$\mathregular{10^{0}}$','$\mathregular{10^{1}}$', '$\mathregular{10^{2}}$', '$\mathregular{10^{3}}$', '$\mathregular{10^{4}}$', '$\mathregular{10^{5}}$'], fontsize = 15)
axarr[0].set_yticklabels(['1','$\mathregular{10^{-1}}$', '$\mathregular{10^{-0.7}}$', '$\mathregular{10^{-0.4}}$', '$\mathregular{10^{-0.2}}$','$\mathregular{10^{-0.09}}$', '$\mathregular{10^{0}}$'], fontsize = 15)
axarr[3].set_yticklabels(['1','$\mathregular{10^{-1}}$', '$\mathregular{10^{-0.7}}$', '$\mathregular{10^{-0.4}}$', '$\mathregular{10^{-0.2}}$','$\mathregular{10^{-0.09}}$', '$\mathregular{10^{0}}$'], fontsize = 15)
axarr[6].set_yticklabels(['1','$\mathregular{10^{-1}}$', '$\mathregular{10^{-0.7}}$', '$\mathregular{10^{-0.4}}$', '$\mathregular{10^{-0.2}}$','$\mathregular{10^{-0.09}}$', '$\mathregular{10^{0}}$'], fontsize = 15)
axarr[0].locator_params(axis='x', nbins=7)


axarr[0].yaxis.set_tick_params(labelsize = 20)
axarr[8].xaxis.set_tick_params(labelsize=20)
axarr[3].yaxis.set_tick_params(labelsize = 20)
axarr[6].xaxis.set_tick_params(labelsize=20)
axarr[6].yaxis.set_tick_params(labelsize = 20)
axarr[7].xaxis.set_tick_params(labelsize=20)
axarr[2].legend(loc = 1, markerscale = 2, prop={'size': 30})


axarr[5].plot(x0[:80], line0[:80], c = 'k', linestyle = '--', linewidth = 3)
axarr[5].plot(x1[:80], line1[:80], c = 'k', linestyle = '--', linewidth = 3)
axarr[5].plot(x2[:130], line2[:130], c = 'k', linestyle = '--', linewidth = 3)


os.chdir(Ffig)
plt.savefig('avalancheduration.png')


    
#============================================================================
#============================================================================
#save video
#---------------
avnum = 9

os.chdir(Ffig)

for i in range(50):
    fig, ax = plt.subplots(figsize= (12,12))
    master = plt.scatter(Fishcoordz[:,0], Fishcoordz[:,1], s=20, c = 'k', alpha = 0.1)
    dotplot = plt.scatter(Fishcoordz[avcoord[i]][:,0], Fishcoordz[avcoord[i]][:,1], s=20, c = 'r', alpha = 1)
    fig.gca().set_aspect('equal', adjustable='box')
    plt.savefig('av' + str(i) + '.tiff')


from skimage import io
from PIL import Image

os.chdir(Ffig)
ogli = sorted(glob.glob('*.tiff'))
omlist = []

for i in range(len(ogli)-1):
    readme = io.imread(Ffig + 'av' + str(i) + '.tiff')
    omlist.append(Image.fromarray(readme))
    #omlist.append(Image.fromarray(rotimglist[i]).convert('L'))
    
omlist[0].save(Ffig + "avmovie.tif", save_all=True,
               append_images=omlist[1:])

#============================================================================
#============================================================================
#PROCESS
#============================================================================
#============================================================================

#PARALLEL PROCESS - without function
#-------------------------------------

#Calculate nearest neighbours for each cell 
#rng = nearest number of cells with which to build neighbour graph from
#dim = define distance of each pixel in x,y,z, cnt = select which % of neighbours to include
#cnt = 0.06, 0.14
#-------------------------------------------------------------------------------------------------
from multiprocessing import Pool
savepath = F10t
processes = 4
pool = Pool(processes)

paramlist = list(range(4))
count = 0
for i in range(np.int(len(coordlist)/processes)):
    paramlist = [Fdrop, experiment, coordlist[count:count+1], 6000, [.8, .8, 15], 0.1], [Fdrop, experiment, coordlist[count+1:count+2], 6000, [.8, .8, 15], 0.1], [Fdrop, experiment, coordlist[count+2:count+3], 6000, [.8, .8, 15], 0.1], [Fdrop, experiment, coordlist[count+3:count+4], 6000, [.8, .8, 15], 0.1]
    pool.starmap(crfn.neighbour, [(paramlist[0]), (paramlist[1]),(paramlist[2]), (paramlist[3])])
    if i == 4:
        paramlist = [Fdrop, experiment, coordlist[20:21], 6000, [.8, .8, 15], 0.1], [Fdrop, experiment, coordlist[21:22], 6000, [.8, .8, 15], 0.14]
        pool.starmap(crfn.neighbour, [(paramlist[0]), (paramlist[1])])
    count+=4
    
#============================================================================
#============================================================================
#Minimum deviation from 2
from sklearn.metrics import mean_squared_error
import heapq

shape = np.load(sparli[0]).shape[0]
expsum = np.zeros((len(sparli),shape))
error = np.full((shape), 1000.0)
extrue = np.full((11), 2.0)
for i in range(shape):
    for t in range(len(sparli)):
        expsum[t,i] = np.load(sparli[t])[i,2]
    if np.isnan(np.sum(expsum)) == False:
        error[i] = mean_squared_error(extrue,expsum[:,i])
        
loc1 = np.where(error == heapq.nsmallest(1, error)[-1]) #smallest
loc2 = np.where(error == heapq.nsmallest(8, error)[-1]) #second smallest

#Minimum variance
from sklearn.metrics import mean_squared_error
import heapq

shape = (20*8)
expsum = np.zeros((len(sparli),shape))
var = np.full((shape), 1000.0)

for i in range(shape):
    for t in range(len(sparli)):
        expsum[t,i] = np.load(sparli[t])[:newshape][i,2]
    if np.isnan(np.sum(expsum)) == False:
        var[i] = np.var(expsum[:,i])
        
loc1 = np.where(var == heapq.nsmallest(1, var)[-1]) #smallest
loc2 = np.where(var == heapq.nsmallest(2, var)[-1]) #smallest
loc3 = np.where(var == heapq.nsmallest(3, var)[-1]) #smallest



#ORDER NAMES
#-----------------
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
    
    
#=============================    
#=============================
#PLOT AVALANCHES    
#=============================
#=============================
# Create datalists
#---------------------------------------------------------------------------
experiment = 'PTZ-WILDTYPE'
num = '02'
os.chdir(Fdrop + 'Project/' + experiment)
binlist = sorted(glob.glob('*E-' + num + '*cutbinarised.npy'))
nnblist = sorted(glob.glob('*E-' + num + '*_nnb.npy')) 
pkglist =  sorted(glob.glob('*E-' + num +'*pkg*')) 
avlist =  sorted(glob.glob('*E-' + num + '*av*'))


#Create time series plot for different avalanches
import scipy.sparse as sp
pkg = p20
pkgs = sp.csr_matrix(p20)
mark = np.unique(pkg, return_counts=True)[0][np.unique(pkg, return_counts=True)[1] > 5][1:] #RETAIN THIS LIST FOR PLOTTING REFERENCE - avalanches with more than 3 cells (can be 3 cells in one time frame, or 3 time frames one cell)
avtime = list(range(mark.shape[0])) #list - each element = 2xarray - 1d: number of cells, 2d: frame number

for t in range(mark.shape[0]):
    avt = np.unique(sp.find(pkgs == mark[t])[1]) #unique time points
    avary = np.zeros((2, len(avt))) #2xarray - 1d: number of cells, 2d: frame number
    count = 0
    for e in avt:
        findmark = pkgs == mark[t]
        avary[0,count] = len(sp.find(findmark)[0][np.where(sp.find(findmark)[1] == e)])
        avary[1,count] = e
        count+=1
    avtime[t] = avary #avtime list is ordered from lowest avalanche number to highest (REMEMBER ONLY TAKING AVALANCHES WITH MORE THAN 2 CELLS)
    
    
#Identify plotting interval
amount = '[:,:3000]'
#find range of avalanches in chosen plotting range
for e in range(len(avtime)):
    if np.max(avtime[e][1]) > 3000:
        cutoff = e
        break 
        
# PLOT - average whole brain fluorescence
#----------------------------------------
experiment = 'PTZ-WILDTYPE'
num = '02'
os.chdir(Fdrop + 'Project/' + experiment)
deltalist = sorted(glob.glob('*E-02*PTZ20*_deltaff.npy'))
f, axarr = plt.subplots(figsize = (20,5))
f.subplots_adjust(hspace=0)
cut = eval('np.load(deltalist[0])' + amount)
maxiar = []
for i in range(cutoff):
    plustime = np.min(avtime[i][1])
    #plt.yscale('log')
    plt.plot(np.append(np.zeros(np.int(plustime)),avtime[i][0]), alpha = 0.5)
    maxiar = np.append(maxiar, avtime[i][0])
    
average = (np.apply_along_axis(np.mean, 0, cut))
normav = average * (np.max(maxiar)/np.max(average))
#plt.yscale('log')
plt.plot(normav)
plt.show()


mini = 500
count = 0
avvec = []
#find range of avalanches in chosen plotting range
for e in range(len(avtime)):
    if np.sum(avtime[e][0]) > mini:
        avvec = np.append(avvec, e)
        count+=1
cutoff = count
print(count)


#PLOT AVALANCHEs
#--------------
avnum = 15
markme = mark[np.int(avvec[avnum])]
time = np.where(pkg == markme)[1]
cells =  np.where(pkg == markme)[0]
coordz = np.load(coordlist[0])
fig, ax = plt.subplots(figsize= (10,10))
master = plt.scatter(coordz[:,0], coordz[:,1], s=35, c = 'k', alpha = 0.1)
dotplot = plt.scatter(coordz[cells][:,0], coordz[cells][:,1], c = time, cmap = 'Spectral_r', s=35, alpha = 1)
fig.colorbar(dotplot, ax = None)
fig.gca().set_aspect('equal', adjustable='box')
os.chdir(Ffig)
plt.savefig('avalanche_large.svg', transparent = True)
plt.show()
