#=======================================================================
def fish_ras(Fishat, fig1, fig2): # Select which fish data to visualise
#=======================================================================
    from matplotlib import pyplot as plt
    import numpy as np
    
# Plotting regime
#--------------------------------------------------------------------------
    f, axarr = plt.subplots(1,2, sharey=True, sharex=True, figsize = (25,25))
    f.subplots_adjust(hspace=0)
    #axarr[0].matshow(Fishat[fig1], cmap = 'Greys', vmin= np.min(Fishat[fig2]), vmax =np.max(Fishat[fig2]))
    axarr[0].set_title(fig1, size = 20)
    axarr[0].set_xlabel('Time', size = 20)
    axarr[0].set_ylabel('Cells', size = 20)
    #axarr[1].matshow(Fishat[fig2], cmap =  vmin= np.min(Fishat[fig2]), vmax =np.max(Fishat[fig2]))
    axarr[1].set_title(fig2, size = 20)
    axarr[1].set_xlabel('Time', size = 20)
    axarr[1].set_ylabel('Cells', size = 20)
    plt.show()

    
#=======================================================================
def fish_zoom(Fishat, figz, a, b, c, d): # Select which fish data to visualise
#=======================================================================
    from matplotlib import pyplot as plt
    import numpy as np
    
    # Zoom into plot
    #---------------------------------------------------------
    plt.figure(figsize=(20,20))
    plt.matshow(Fishat[figz][a:b,c:d], fignum=1)
    plt.title(figz, size = 20)
    plt.ylabel("Cells", size = 15)
    plt.xlabel("Time", size = 15)
    plt.show()
    
    
    
    
#=======================================================================
def fish_mamp(Fishat, fli, okgo, bins, meantime, meancells, plt1, plt2, plt3): # Select which fish data to visualise
#=======================================================================

    import numpy as np
    from matplotlib import pyplot as plt
    import pandas as pd
    import ptitprince as pt
    import seaborn as sns
    import os    
    
    okgoyo = list(range(len(okgo)))
    
    for i in range(len(okgo)):
        okgoyo[i] = Fishat[okgo[i]]
    
    
    if meantime == 'yes':
    
    #Whole brain - single cell, all averaged across recording in specific time bins 
    #Single fish - time binned amplitude over recording, with all cells 
    #averaged together at each time point
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------

        newvlist = list(range(len(okgoyo)))    
        for y in range(len(okgoyo)):
            newv = np.zeros(shape = (okgoyo[y].shape[0],int(okgoyo[y].shape[1]/bins)))
            print("Binning fish " + str(y))
            for i in range(okgoyo[y].shape[0]):
                
                cells = okgoyo[y][i,:]
                reshape = np.reshape(cells, (int(len(cells)/bins), bins))
                meanamp = np.apply_along_axis(np.mean, 1, reshape)
                newv[i] = meanamp
                mymean = np.apply_along_axis(np.mean, 0, newv)
                newvlist[y] = mymean
    
    
    
    #Plotting regime for amplitude trajectory 
    #-------------------------------------------------------------------------------        
        plt.figure(figsize=(10,3))
        plt.plot(newvlist[plt1], label = 'baseline')
        plt.plot(newvlist[plt2], label = 'ptz 5mm')
        plt.plot(newvlist[plt3], label = 'ptz 20mm')
        plt.ylabel("Mean amplitude (all cells)", size = 10)
        plt.xlabel("Time bins", size = 10)
        plt.title('Amplitude trajectory of pre/ictal states', size = 15)
        plt.legend(loc='upper left')
        plt.savefig('trajectory')
        plt.show()
            
    #Create pandas dataframe from mean bin amplitude and plot as violin plot
    #-------------------------------------------------------------------------------
        dit = {'baseline': newvlist[0], 'P5': newvlist[1], 'P20': newvlist[2]}
        ditf = pd.DataFrame.from_dict(dit, orient = 'index')
        trudit = ditf.transpose()
        sns.set()
        
        fig2 = plt.figure(figsize=(7,7))
        seal2 = sns.violinplot(data=trudit)
        seal2.set_title('Distribution of amplitudes across time bins', size = 15)
        seal2.set_ylabel("Mean amplitude (all cells)", size = 15)
        plt.show()
    
        Fishab = {}
    
        for i in range(len(okgo)):
            Fishab.update({okgo[i]: newvlist[i]})
    
    
    if meancells == 'yes':
        
    # Compare mean amplitude for baseline, 5mm and 20mm ptz (one fish), 
    # using newvlist mean data from above plot
    #-----------------------------------------------------------------------
    #-----------------------------------------------------------------------
    
        #sns.set(style="darkgrid")
        #sns.set(style="whitegrid")
        #sns.set_style("white")
        #sns.set(style="whitegrid",font_scale=2)
        
        patalist = list(range(len(okgoyo)))
        filtlist = list(range(len(okgoyo)))
        afc = list(range(len(okgoyo)))

        # Calculate mean amplitude of each cells per condition 
        #---------------------------------------------------
        for t in range(len(okgoyo)):
            meanz = list(range(okgoyo[t].shape[0]))
            for j in range(okgoyo[t].shape[0]):
                ofc = np.mean(okgoyo[t][j])
                meanz[j] = ofc
            afc[t] = meanz 

        dat = {'baseline': afc[0], 'P5': afc[1], 'P20': afc[2]}
        datf = pd.DataFrame.from_dict(dat, orient='index')
        trudat = datf.transpose()
        sns.set()

        #Raincloud plot
        #------------------------
        pal = "Set2"; sigma = .2
        ax=pt.RainCloud(data = trudat, palette = pal, bw = sigma,
                 width_viol = .8, figsize = (15,10), move = 0.2, alpha = 1, orient = 'h')
        ax.set_title('Mean amplitude per cell', size = 10)
        plt.show()
        
        #Violin plot showing distribution of amplitude for each cell 
        #-------------------------------------------------------------------------------
        fig1 = plt.figure(figsize=(7,7))
        seal1 = sns.violinplot(data=trudat)
        seal1.set_ylabel("Mean amplitude per cell", size = 15)
        seal1.set_title('Distribution of cell amplitudes of pre/ictal states', size = 15)
        plt.show()

        
        Fishacl = {}
        for i in range(len(okgo)):
            Fishacl.update({okgo[i]: afc[i]})
    
    if meancells == 'yes' and meantime == 'yes':
        return(Fishab, Fishacl)
    
    if meancells == 'no':
        return(Fishab, {})
    if meantime == 'no':
        return({},Fishacl)
        
    
#=======================================================================
def fish_peaks(okgo, fli, Fishat, std, py, smooth, highcut, lowcut, stdlim, exf, exc): # Select which fish data to visualise
#=======================================================================
    import numpy as np
    from scipy.fftpack import rfft, irfft, fftfreq
    import matplotlib.pyplot as plt
    from scipy.signal import find_peaks
    import numpy as np
    from scipy import stats
    from scipy import signal
    
    inputlist = list(range(len(okgo)))
    oinputlist = list(range(len(okgo)))
    pkblist = list(range(len(okgo)))
    pkmamp = list(range(len(okgo)))
    pypeaklist = list(range(len(okgo)))
    wherepeaklist = list(range(len(okgo)))
    peakamp = list(range(len(okgo)))
    
    #Detect peaks - smooth then peak detect
    #--------------------------------------------  

    okgoyo = list(range(len(okgo)))    
    for i in range(len(okgo)):
        okgoyo[i] = Fishat[okgo[i]]
    
    inputlist = list(range(len(okgo)))
    
    
    if smooth == 'yes':
    
    # Smooth traces - low pass filter
    #-----------------------------------------
        for y in range(len(okgoyo)):
            filterme = np.zeros(okgoyo[y].shape)
    
            # loop through each time point and apply filter
            #------------------------------------------
            for i in range(okgoyo[y].shape[0]):
                d = okgoyo[y][i,:]
                f_signal = rfft(d)
                f_signal[0:(2*lowcut+1)] = 0
                f_signal[(2*highcut+1):len(f_signal)] = 0
                filterme[i,:] = irfft(f_signal)
            
            inputlist[y] = filterme
        
    if smooth == 'no':
        
        for i in range(len(okgoyo)):
            inputlist[i] = okgoyo[i]
    
    if py == 'yes':

    #Detect peaks
    #--------------------------------------------------------------------------   
        pylist = list(range(len(okgoyo)))

        for y in range(len(okgoyo)):
    
            peakz = list(range(inputlist[y].shape[0]))
            pks = np.zeros(okgoyo[y].shape)
            pokeamp = list(range(inputlist[y].shape[0]))
            
            for i in range(inputlist[y].shape[0]):
                peaks, _ = find_peaks(inputlist[y][i], height=0)
                peakz[i] = peaks
                pokeamp[i] = okgoyo[y][i,peaks]
            pypeaklist[y] = peakz
            pkblist[y] = pks
            peakamp[y] = pokeamp
         
        Fishpyp = {}
        for y in range(len(okgo)):
            Fishpyp.update({okgo[y]: peakamp[y]})   #dict for peak amplitude values for each cell
      
        
        #Show peaks
        #--------------------------------------------------------------------------   
        x = inputlist[exf][exc]
        peaky = pypeaklist[exf][exc]
        plt.plot(x)
        plt.plot(peaky, x[peaky], "x")
        plt.plot(np.zeros_like(x), "--", color="gray")
        plt.show()
        
        
    if std == 'yes':
        
        for i in range(len(okgoyo)):
            oinputlist[i] = okgoyo[i]
        for y in range(len(okgoyo)):

        # Find activity peaks
        #---------------------------------------------------------------------------
            pks = np.zeros(okgoyo[y].shape)
            pkampl = list(range(okgoyo[y].shape[0]))
            wherepeaks = list(range(okgoyo[y].shape[0]))
        
            for i in range(okgoyo[y].shape[0]):
                d = oinputlist[y][i,:]                                                            
                sem = np.std(d)
                p   = np.where(d > stdlim*sem)[0]   #p is indeces of all peaks
                pks[i,p] = 1
                pkampl[i] = okgoyo[y][i,p]
                wherepeaks[i] = p
    
            pkblist[y] = pks    #binary list of where peaks are
            pkmamp[y] = pkampl  #peak amplitude values
            wherepeaklist[y] = wherepeaks
       
        Fishpkb = {}
        Fishstdp = {}
        for i in range(len(okgo)):
            Fishpkb.update({okgo[i]: pkblist[i]})
            Fishstdp.update({okgo[i]: pkmamp[i]})
        
        
        
        #Show peaks
        #--------------------------------------------------------------------------   
        x = inputlist[exf][exc]
        peaky = wherepeaklist[exf][exc]
        plt.plot(x)
        plt.plot(peaky, x[peaky], "x")
        plt.plot(np.zeros_like(x), "--", color="gray")
        plt.show()

    
    if py == 'yes' and std == 'yes':
        return(Fishpyp, Fishpkb, Fishstdp)
    
    if py == 'no' and std == 'yes':
        return({}, Fishpkb, Fishstdp)
    
    if py == 'yes' and std == 'no':
        return(Fishpyp, {},{})
    

#=======================================================================
def fish_pkamps(Ffig, Fishat, Fishpyp, Fishstdp, fli, okgo, typ, rain , viol, plt1, plt2, plt3): # Select which fish data to visualise
#=======================================================================
    
    import numpy as np
    import os
    from scipy.fftpack import rfft, irfft, fftfreq
    import matplotlib.pyplot as plt
    from scipy.signal import find_peaks
    import numpy as np
    from scipy import stats
    from scipy import signal
    import pandas as pd
    import ptitprince as pt
    import seaborn as sns
    
    afcp = list(range(len(okgo)))

    #Calculate peak mean amplitude for each cell in each condition
    #--------------------------------------------------------------------
    for t in range(len(okgo)):
        

        pkmeanz = list(range(Fishat[okgo[t]].shape[0]))
        
        if typ == 'si':
            okpoyo = list(range(len(okgo))) 
            for i in range(len(okgo)):
                okpoyo[i] = Fishpyp[okgo[i]]
            var = okpoyo   
    
        elif typ == 'mu':
            oksoyo = list(range(len(okgo)))
            for i in range(len(okgo)):
                oksoyo[i] = Fishstdp[okgo[i]]
            var = oksoyo
          
              
        for j in range(Fishat[okgo[t]].shape[0]):
            ofcp = []
            ofcp = np.mean(var[t][j])
            pkmeanz[j] = ofcp
    
        afcp[t] = pkmeanz
    
    Fish1cpm = {}
    for i in range(len(okgo)):
        Fish1cpm.update({okgo[i]: afcp[i]})
        
    
   
    #Plot 
    #--------------
    dut = {'baseline': afcp[plt1], 'PTZ 5mM': afcp[plt2], 'PTZ 20mM': afcp[plt3]}
    dutf = pd.DataFrame.from_dict(dut, orient='index')
    trudut = dutf.transpose()
    sns.set()

    if rain == 'yes':
    #Raincloud plot
    #------------------------
        pal = "Set2"; sigma = .2
        ax=pt.RainCloud(data = trudut, palette = pal, bw = sigma, width_viol = .8, figsize = (15,10), move = 0.2, alpha = 1, orient = 'h')

        ax.set_title('Mean amplitude per cell', size = 10)
        plt.savefig('Rainpeakamp')
        plt.show()
        

    if viol == 'yes':
              
        plt.figure(figsize=(10,10))
        ax = sns.violinplot(data=trudut, palette=['g','b','r'])
        ax.set_ylabel("Mean amplitude [au]", size = 20, color = 'black')
        plt.xticks(np.arange(3), ('Baseline','PTZ 5mM', 'PTZ 20 mM'), size = 15, color = 'black')
        plt.yticks(color = 'black')
        #ax.set_title('Single cell firing event amplitude', size = 15)
        os.chdir(Ffig)
        plt.savefig('Violpeakamp.png')
        plt.show()
        
              
    return(Fish1cpm)


