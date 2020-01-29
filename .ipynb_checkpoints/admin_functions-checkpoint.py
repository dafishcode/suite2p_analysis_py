#SORT
#=============================
#=============================

#=============================
def name(path, experiment, start, end): #add 0 to number if above 10
#=============================
    import os 
    count = 0
    listme = list(range(start, end+1))
    os.chdir(path + 'Project/' + experiment)
    for i in range(start, end+1):
        if i < 10:
            num = '0' + str(i)
        elif i >9:
            num = str(i)
        listme[count] = num
        count+=1
    return(listme)

#=============================
def nameli(path, experiment, inp, string): #return name list
#=============================
    import os 
    import glob
    os.chdir(path + 'Project/' + experiment)
    if inp < 10:
        out = '0' + str(inp)
    elif inp >9:
        out = str(inp)
    return(sorted(glob.glob('*E-' + str(out) + string)))


#==============================
def save_name(i, name_li): #find save name
#===============================
    return(name_li[i][:name_li[i].find('run')+6])

#PROCESS
#=============================
#=============================

#=====================================================================
def parallel(cores, datalist1, datalist2, func, paramlist): #make sure n_cores is divisible by total number
#=====================================================================
    from multiprocessing import Pool
    import numpy as np
    pool = Pool(cores)
    paramlist_levels = list(range(cores))
    count = 0
    for i in range((np.int(len(datalist1)/cores))):
        for e in range(len(paramlist_levels)):
            paramlist_levels[e] = np.append(paramlist, datalist1[count], datalist2[count])
            count+=1
        output = pool.starmap(func, paramlist_levels)
    
#MATHS
#=============================
#=============================

