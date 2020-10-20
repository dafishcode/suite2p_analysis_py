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
        kmeans = KMeans(n_clusters=16, random_state=0).fit(corr)
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


#=====================
#=====================
class ws_netsim: 
#=====================
#=====================
    """
    Class to build watts-strogatz networks and run avalanche simulations
    dist = distance matrix between all nodes in network
    
    """

    #========================
    def __init__(self,dist):
    #========================
        import numpy as np
        self.dist = dist
    

    #BUILD NETWORK
    #=================
    #=================
    #====================================
    def k_neighbours(self, edge_density, mode):
    #====================================
        import numpy as np
        """
        Form connections with k-nearest neighbours
        
        edge_density = number of neighbours to connect to
        mode = directed or undirected
        """
        
        # Loop through rows of distance matrix to find k_neighbours
        #-----------------------------------------------------------------------------
        for row in range(self.dist.shape[0]):
            k_neighbours = int(edge_density) #Find k_neighbours for each cell
            neighbours = self.dist[row,].argsort()[:k_neighbours+1][::-1] #find neighbours 
            self.A[row,neighbours[:-1]] = 1 #make all edges of neighbours connected in network
            
            if mode == 'undirected':
                self.A[neighbours[:-1],row] = 1
        return(self)

    #=====================================
    def net_generate(self, edge_density, p, mode):
    #=====================================
        """
        Generate random small world graph with specific Edge density. The Watts-Strogatz model has (i) a small average shortest path length, and (ii) a large clustering coefficient. The algorithm works by assigning a pre-defined number of connections between k-nearest neighbours - it then loops through each node and according to some uniform probability re-assigns its edges from its connected k-nearest neighbours and a random unconnected node. 

            edge_density = number of k_nearest neighbours each node is connected to
            p = probability of an edge being randomly re-assigned
            mode = directed or undirected
        """
        import numpy as np
        import networkx as nx
        import random
        import copy
        
        if mode!= 'directed' and mode!= 'undirected': 
            print('Select directed or undirected')
            exit()
        self.A = np.zeros(self.dist.shape)
        self.k_neighbours(edge_density, mode)

        # Rewire connections with certain probability
        #-----------------------------------------------------------------------------
        
        if mode == 'undirected':
            [rows, cols]    = np.where(np.triu(self.A) == 1) 
            probs           = np.random.uniform(size = rows.shape[0]) #Generate random values for each connection 
            edges_to_change = np.where(probs <= p)[0] #see which values are randomly changed
            
            for e in range(edges_to_change.shape[0]): #Loop through edges to change
                this_edge = edges_to_change[e]
                self.A[rows[this_edge], cols[this_edge]] = 0         # switch off old edge
                self.A[cols[this_edge], rows[this_edge]] = 0

                where_0 = np.where(self.A[rows[this_edge]] == 0)[0] #find possible connections to reassign to
                new_edge = random.choice(where_0[np.where(where_0 !=rows[this_edge])[0]]) #randomly choose one - ignoring any connections on the diagonal 
                #Assign connection
                self.A[rows[this_edge], new_edge] = 1        # switch on new edge
                self.A[new_edge, rows[this_edge]] = 1
        
        if mode == 'directed':
            [rows, cols]    = np.where(self.A == 1) 
            probs           = np.random.uniform(size = rows.shape[0]) #Generate random values for each connection 
            edges_to_change = np.where(probs <= p)[0] #see which values are randomly changed
        
            # Rewire connections with certain probability
            #-----------------------------------------------------------------------------
            [rows, cols]    = np.where(self.A == 1) 
            probs           = np.random.uniform(size = rows.shape[0]) #Generate random values for each connection 
            edges_to_change = np.where(probs <= p)[0] #see which values are randomly changed

            for e in range(edges_to_change.shape[0]): #Loop through edges to change
                this_edge = edges_to_change[e]
                self.A[rows[this_edge], cols[this_edge]] = 0         # switch off old edge

                where_0 = np.where(self.A[rows[this_edge]] == 0)[0] #find possible connections to reassign to
                new_edge = random.choice(where_0[np.where(where_0 !=rows[this_edge])[0]]) #randomly choose one - ignoring any connections on the diagonal 
                #Assign connection
                self.A[rows[this_edge], new_edge] = 1        # switch on new edge
        return(self)

    
    #CALCULATE CYCLES
    #=================
    #=================
    #===========================
    def cycles_calculate(self, edge_density, p, mode):
    #===========================
        import networkx as nx
        import numpy as np
        
        cyc_mat = self.net_generate(edge_density, p, mode).A #matrix to calculate cycles
        G = nx.from_numpy_matrix(cyc_mat)
        cyc = nx.algorithms.cycle_basis(G)
        edge =  int(np.sum(cyc_mat))
        self.cycles = len(cyc)
        self.edges = edge
        return(self)
        
    #===========================
    def cycles_median(self, edge_density, p, n_samp, mode):
    #===========================
    #select median cycles number for simulations - ensure you capture non-skewed cycle values
        import networkx as nx
        import numpy as np
        cyc_list = list(range(n_samp)) #list containing cycle densities for each iteration
        cyc_mat_list = list(range(n_samp)) #list containing each generated matrix
        for i in range(n_samp):
            curr_mat = self.net_generate(edge_density, p, mode).A #matrix to calculate cycles
            G = nx.from_numpy_matrix(curr_mat)
            cyc = nx.algorithms.cycle_basis(G)
            edge =  int(np.sum(curr_mat))
            cyc_mat_list[i] = curr_mat
            cyc_list[i] = len(cyc)/edge
        if n_samp % 2 == 0:
            self.sim_A  = cyc_mat_list[min(range(len(cyc_list)), key=lambda x: abs(cyc_list[x]-np.median(cyc_list)))] #matrix to run simulation on

        else:
            self.sim_A  = cyc_mat_list[np.where(cyc_list == np.median(cyc_list))[0][0]] #matrix to run simulation on

        return(self) 
    
    
    #BUILD WEIGHT MATRIX
    #===================
    #===================
    # Simple sigmoid function to 'soften' the exponential
    #===========================
    def sig(self, x):
    #===========================
        import numpy as np
        self.sig_output = 1 / (1+np.exp(-x))
        return(self)
    
    # Conversion from distance to edge weights, scaled (itself exponentially) by s
    #====================================
    def dist2edge(self, distance, divisor, soften, s):
    #===================================
        import numpy as np
        self.edge_weight_out = np.exp(s/5)*self.sig(np.exp(-soften/np.exp(s)*distance)).sig_output/divisor
        return(self)  
    
    #===========================
    def adjmat_generate(self, s, edge_density, p, n_samp, divisor, soften, mode):
    #===========================
        import numpy as np
        import copy
        mat = np.zeros((self.dist.shape))
        
        curr_mat = self.cycles_median(edge_density, p, n_samp, mode).sim_A
        
        [rows, cols]    = np.where(np.triu(curr_mat) == 1) 
        for e in range(len(rows)):
            edge_weight = self.dist2edge(self.dist[rows[e], cols[e]], divisor, soften, s).edge_weight_out
            mat[rows[e], cols[e]] = edge_weight 
            mat[cols[e], rows[e]] = edge_weight
        self.adj_mat = copy.deepcopy(mat)
            
        return(self)
    
    
    
    #SIMULATE AVALANCHES
    #===================
    #===================
    
    #Find cells to propagate
    #=====================================================
    def propagate_neighbours(self, curr_mat, start_node):
    #=====================================================
        import numpy as np
        self.prop_nodes = []
        nodes = np.where(curr_mat[start_node] > 0) [0]
        weights = curr_mat[start_node][nodes]
        for f in range(len(nodes)):
            if weights[f] > np.random.uniform(0, 1):
                self.prop_nodes = np.append(self.prop_nodes, nodes[f])
        return(self)

    
    #Simulate 
    #===========================
    def simulate(self,  s, edge_density, p, n_samp, divisor, soften, cutoff, n_sims, mode):
    #===========================
        import numpy as np
        curr_mat = self.adjmat_generate(s, edge_density, p, n_samp, divisor, soften, mode).adj_mat

        self.av_size = []
        self.av_dur = []

        for i in range(n_sims):
            #Decide start node
            start_node = np.random.uniform(0, curr_mat.shape[0]-1)
            down = int(start_node)
            up= int(start_node)+1
            if np.random.uniform(down, up) >= start_node:
                start_node = up
            else:
                start_node = down


            #Initialise avalanche - ping first node
            t_nodes = self.propagate_neighbours(curr_mat, start_node).prop_nodes #Find connected neighbours > threshold
            curr_list = t_nodes
            iterate = 'yes'

            if len(t_nodes) > 1: #must have at least 3 cells to begin avalanche
                all_nodes = np.append(start_node, t_nodes)
                timesteps = 1

                while iterate == 'yes':
                    tplus_nodes = []
                    for z in range(len(curr_list)):
                        #List of all nodes active in next timestep
                        tplus_nodes = np.append(tplus_nodes, self.propagate_neighbours(curr_mat, int(curr_list[z])).prop_nodes)

                    all_nodes = np.append(all_nodes, tplus_nodes)
                    timesteps+=1
                    curr_list = tplus_nodes

                    if len(all_nodes) > cutoff:
                        iterate = 'no'

                    if len(tplus_nodes) == 0: #if no more active cells - stop
                        iterate = 'no'


                self.av_size = np.append(self.av_size, len(all_nodes)) 
                self.av_dur = np.append(self.av_dur, timesteps)

            else:
                continue

        return(self)
    
    
    
    
    
#=====================
#=====================
class ba_netsim: 
#=====================
#=====================
    """
    Class to build barabasi-albert networks and run avalanche simulations
    dist = distance matrix between all nodes in network
    """

    #========================
    def __init__(self,dist):
    #========================
        import numpy as np
        self.dist = dist
    

    #BUILD NETWORK
    #=================
    #=================
    
    #=====================================
    def sample(self, seq, m):
    #=====================================
        """ Return m unique elements from seq.
        """
        import random
        import numpy as np
        
        #make targets a set - only contains unique elements
        targets=set()
        while len(targets)<m:
            x=random.choice(seq)
            targets.add(x) #add method only adds if x is not already in target set
        return np.array(list(targets))

    
    #=====================================
    def net_generate(self, edge_density):
    #=====================================
        """
        Generate Barabasi-Albert preferential attachment network. BA model starts with k initial nodes, and k edges 
        - each new node will connect to k nodes with p(number of edges already connected to each node). 
        
            edge_density = number of edges of each node
            
        """
        
        import numpy as np
        import networkx as nx
        import random
        import copy
        
        self.A = np.zeros(self.dist.shape) #initialise graph

        current_n = edge_density #current number of nodes

        # Nodes to connect to from current node
        nodes_out =list(range(edge_density))

        # Sequence of all nodes connected (in and out) - can sample from this 
        node_counts=[]

        #iterate until number of nodes = n
        while current_n < self.dist.shape[0]:
            self.A[current_n, nodes_out] = 1
            self.A[nodes_out, current_n] = 1 

            #add current nodes receiving outgoing connections to node sequence
            node_counts.extend(nodes_out)

            #list of incoming connections for current node - i.e. repeated sequence of current node
            nodes_in = [current_n]*edge_density

            #add current node (as many times as it sends out connections - assumes undirected network) to node sequence
            node_counts.extend(nodes_in)

            #update nodes_out - uniformly sample from sequence of node_counts
            nodes_out = self.sample(node_counts, edge_density)

            current_n +=1
        return(self)
    
    
    
    #CALCULATE CYCLES
    #=================
    #=================
    #===========================
    def cycles_calculate(self, edge_density):
    #===========================
        import networkx as nx
        import numpy as np
        
        cyc_mat = self.net_generate(edge_density).A #matrix to calculate cycles
        G = nx.from_numpy_matrix(cyc_mat)
        cyc = nx.algorithms.cycle_basis(G)
        edge =  int(np.sum(cyc_mat))
        self.cycles = len(cyc)
        self.edges = edge
        return(self)
    
    
    #BUILD WEIGHT MATRIX
    #===================
    #===================
    # Simple sigmoid function to 'soften' the exponential
    #===========================
    def sig(self, x):
    #===========================
        import numpy as np
        self.sig_output = 1 / (1+np.exp(-x))
        return(self)
    
    # Conversion from distance to edge weights, scaled (itself exponentially) by s
    #====================================
    def dist2edge(self, distance, divisor, soften, s):
    #===================================
        import numpy as np
        self.edge_weight_out = np.exp(s/5)*self.sig(np.exp(-soften/np.exp(s)*distance)).sig_output/divisor
        return(self)  
    
    #===========================
    def adjmat_generate(self, s, edge_density, divisor, soften):
    #===========================
        import numpy as np
        import copy
        mat = np.zeros((self.dist.shape))
        
        curr_mat = self.net_generate(edge_density).A #matrix to calculate cycles
        
        [rows, cols]    = np.where(np.triu(curr_mat) == 1) 
        for e in range(len(rows)):
            edge_weight = self.dist2edge(self.dist[rows[e], cols[e]], divisor, soften, s).edge_weight_out
            mat[rows[e], cols[e]] = edge_weight 
            mat[cols[e], rows[e]] = edge_weight
        self.adj_mat = copy.deepcopy(mat)
            
        return(self)
    
    
    
    #SIMULATE AVALANCHES
    #===================
    #===================
    
    #Find cells to propagate
    #=====================================================
    def propagate_neighbours(self, curr_mat, start_node, r_e):
    #=====================================================
        import numpy as np
        self.prop_nodes = []
        nodes = np.where(curr_mat[start_node] > 0) [0]
        weights = curr_mat[start_node][nodes]
        for f in range(len(nodes)):
            if weights[f] > np.random.uniform(0, 1):
                self.prop_nodes = np.append(self.prop_nodes, nodes[f])
                
        if r_e > np.random.uniform(0,1):
            self.prop_nodes = np.append(self.prop_nodes, start_node)
        return(self)

    
    #Simulate 
    #===========================
    def simulate(self,  s, edge_density, max_e, divisor, soften, cutoff, n_sims, iterate):
    #===========================
        import numpy as np
        curr_mat = self.adjmat_generate(s, edge_density, divisor, soften).adj_mat
        degree = sum(self.A)
        degree_scaled = (max_e/max(degree)) * degree 
        
        


        
        self.av_size = []
        self.av_dur = []
        
        #iterate simulation calculation for less-noisy distribution
        for x in range(iterate):
            
            for i in range(n_sims):
                #Decide start node
                start_node = np.random.uniform(0, curr_mat.shape[0]-1)
                down = int(start_node)
                up= int(start_node)+1
                if np.random.uniform(down, up) >= start_node:
                    start_node = up
                else:
                    start_node = down


                #Initialise avalanche - ping first node
                t_nodes = self.propagate_neighbours(curr_mat, start_node, degree_scaled[start_node]).prop_nodes #Find connected neighbours > threshold
                curr_list = t_nodes
                iterate = 'yes'

                if len(t_nodes) > 1: #must have at least 3 cells to begin avalanche
                    all_nodes = np.append(start_node, t_nodes)
                    timesteps = 1

                    while iterate == 'yes':
                        tplus_nodes = []
                        for z in range(len(curr_list)):
                            #List of all nodes active in next timestep
                            tplus_nodes = np.append(tplus_nodes, self.propagate_neighbours(curr_mat, int(curr_list[z]), degree_scaled[int(curr_list[z])]).prop_nodes)

                        all_nodes = np.append(all_nodes, tplus_nodes)
                        timesteps+=1
                        curr_list = tplus_nodes

                        if len(all_nodes) > cutoff:
                            iterate = 'no'

                        if len(tplus_nodes) == 0: #if no more active cells - stop
                            iterate = 'no'


                    self.av_size = np.append(self.av_size, len(all_nodes)) 
                    self.av_dur = np.append(self.av_dur, timesteps)

                else:
                    continue

        return(self)
    