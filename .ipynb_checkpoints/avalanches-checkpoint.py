#PROCESS
#------------
#------------

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
        np.save(Fdrop + 'Project/' + experiment +os.sep + coordlist[i][:coordlist[i].find('run')+6] + '_' + str(cnt) + 'nnb.npy', nnb)
    return(nnblist)


#ANALYSIS
#------------
#------------
#=======================================================================
def avalanche(nnb, bind, savepath, experiment): # duration = yes convergence (no back propagation, earliest avalanche consumes meeting avalanche, and later avalanche terminates), cells in t must be active in t+1)
#=======================================================================
    import numpy as np
    import os
    import itertools

#Calculate avalanche size + duration
#-----------------------------------
    avlist = list(range(len(nnblist)))
    pkglist = list(range(len(nnblist)))
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
        if i% round(10*binarray.shape[1]/100) == 0: print('doing time step ' + str(i) + 'of' + str(binarray.shape[1]) + 'for fish ' + str(y))
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
                
                elif t>30:
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
    
    np.save(savepath + 'Project/' + experiment + os.sep + binlist[y][:binlist[y].find('run')+7] + 'av.npy', av)
    np.save(savepath + 'Project/' + experiment + os.sep + binlist[y][:binlist[y].find('run')+7] + 'pkg.npy', pkg)
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
        fit = powerlaw.Fit(data, xmax = 100)
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
def cellplot(Ftm, Fdrop,F10t, experiment, fnum, prefix, condition, plane, cell, xshift, yshift): # Plot cells and neighbours over image 
#=======================================================================
    import os
    import glob
    import numpy as np
    from matplotlib import pyplot as plt
    cs = []
    
    Freg = Ftm + 'Project/' + experiment + '-' + fnum + prefix
    os.chdir(Freg)
    opslist = sorted(glob.glob('*' + condition + '*' +'plane' + str(plane) + '*ops.npy'))
    ci   = np.where(cs[:,2] == plane)[0]    # Index of plane coordinates in long list
    os.chdir(F10t + 'Project/' + experiment)
    nnblist = sorted(glob.glob('*' + fnum +  '*' + condition + '*0.06nnb.npy'))
    nnb = np.load(nnblist[0])
    os.chdir(Fdrop + 'Project/' + experiment)
    coordlist = sorted(glob.glob('*' + fnum +  '*' + condition + '*realcoord.npy'))
    cs = np.load(coordlist[0])             # 3D array of xyz coordinates



    
    if len(opslist) + len(coordlist) + len(nnblist) >3:
        print('More than one fish image loaded')
        
    # Plot
    #------------------------------------------------------
    ops = np.load(Freg + opslist[0])
    ops = ops[()]
    raw = ops['meanImg']

    # Pull out data from fish structure
    #----------------------------------------------------------------------------------------

    # Actual plotting routines
    #----------------------------------------------------------------------------------------
    plt.figure(figsize = (15,15))
    plt.imshow(raw)
    plt.scatter(cs[ci,0]+xshift, cs[ci,1]-yshift, s = 10, c = nnb[ci[cell],ci], cmap = 'prism')
    plt.show()

            

#================================================
def powerfit_param(Fdrop, F10t, experiment, cutoff, num):
#=================================================
    import os
    import numpy as np
    import glob
    import powerlaw
    
    os.chdir(F10t + 'Project/' + experiment)
    itav = sorted(glob.glob('*E-' + str(num) + '*nnbav.npy*'))
    os.chdir(Fdrop + 'Project/' + experiment)
    coord = sorted(glob.glob('*E-' + str(num) + '*BLN*realcoord*')) 
    cells = np.load(coord[0]).shape[0]
    os.chdir(F10t + 'Project/' + experiment)
    
    paramar = np.zeros((len(itav), 6))
    for i in range(len(itav)):
        data = np.load(itav[i])[1]
        if len(data) < 10 : #if less than 10 avalanches skip 
            continue
        maxi = np.max(np.unique(data, return_counts = True)[0][np.unique(data, return_counts = True)[1] > cutoff]) #fit power law max - maximum value that appears more than 3 times
        fit = powerlaw.Fit(data, discrete = True, xmax = maxi) #fit power law to data - MLE    
        alpha = fit.power_law.alpha
        sigma = fit.power_law.sigma
        if maxi - fit.xmin > 3:
            R, p = fit.distribution_compare('truncated_power_law', 'lognormal', normalized_ratio=True)
            paramar[i,0] = R
            paramar[i,1] = p
        paramar[i,2] = alpha
        paramar[i,3] = sigma
        paramar[i,4] = maxi
        paramar[i,5] = maxi/cells

    np.save(Fdrop + 'Project/' + experiment + os.sep + coord[0][:coord[0].find('run')+6] + '_' + 'durparamsweep1.npy', paramar)
    print('Done fish num ' + str(f))     

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