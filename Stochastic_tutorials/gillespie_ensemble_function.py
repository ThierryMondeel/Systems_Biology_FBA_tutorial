
# coding: utf-8

# In[72]:

def gillespie_ensemble(N_CELLS = 5000, k_burst = 50.0, k_deg = 1.0, burst_size = 1.0 ):
    
    # IMPORT stuff
    import random # 'random' is a module to generate random variables
    import matplotlib.pyplot as plt # package to make plot
    get_ipython().magic(u'matplotlib inline')
    import numpy as np 
    import pickle
    import math 
    import scipy.stats.stats as st
    import pandas as pd
    import operator
    
    # DEFINE A CLASS "CELL" WITH THE FOLLOWING FEATURES
    class CELL:

        def __init__(self, cell_type, k_deg, k_burst, deg_size, burst_size, min_time, initial_mrna):

            self.cell_type = cell_type

            self.k_deg= k_deg
            self.k_burst = k_burst
            self.deg_size = deg_size
            self.burst_size = burst_size

            self.trajectory = {min_time:initial_mrna}

            self.van_Kampen_Fano_factor = ((self.deg_size + self.burst_size)/((self.deg_size)*2+0.00))
            self.macroscopic_mean = ((self.k_burst*self.burst_size)/((self.k_deg*self.deg_size)+0.00))
            # NB I changed the macroscopic mean adding the growth rate to the denominator
            tau = 1/((self.k_deg*self.deg_size)+0.00)
            self.Half_life = round(tau*(np.log(2)),2)

        def features(self):
            return "Cell type {} * k deg: {} * k burst: {} * deg size: {} * burst size: {}".format(self.cell_type, 
                                                                                                      self.k_deg, 
                                                                                                      self.k_burst, 
                                                                                                      self.deg_size, 
                                                                                                      self.burst_size )
    _break_ = "\033[0;34m--------------------------------------------------------------------------------------------------------------- \033[0m"
    _star_break_ = '     **     '
    
        
    # ---------------------------------------------------------- #
       
    N_CELLS = N_CELLS

    Dict_min_time_Fano_var_mean = {}

    min_time = 0

    Distribution_dict_list = {}

    last_time_point = []  
    cell_status = list(range(0,N_CELLS))
    cell_track = list(range(0,N_CELLS))
    # ---------------------------------------------------------- #
    
    
    # PARAMETERS CELL TYPE 

    cell_type = "X"

    k_deg = k_deg

    k_burst = k_burst

    deg_size,burst_size = 1.0, burst_size


    min_time_zero = 0
    initial_mrna_zero = 0
    # ---------------------------------------------------------- #


    for n in range(N_CELLS):

        cell_track[n] =  CELL(cell_type, k_deg, k_burst, deg_size, burst_size, min_time_zero, initial_mrna_zero)
        random_m = (np.random.uniform(0.8,1.2))*cell_track[n].macroscopic_mean
        random_std = (np.random.uniform(0.8,1.2))*round(np.sqrt(cell_track[n].macroscopic_mean*cell_track[n].van_Kampen_Fano_factor),2)
        p = max([round(int(np.random.normal(loc= random_m, scale=random_std))),0])
    #     p = max([round(int(np.random.normal(loc= cell_track[n].macroscopic_mean, scale=round(np.sqrt(cell_track[n].macroscopic_mean*cell_track[n].van_Kampen_Fano_factor),2)))),0])


        cell_track[n].trajectory = {0: p} # {time zero: poisson mRNA}
        cell_status[n] = p
        last_time_point.extend([0])


        Distribution_dict_list[0] = cell_status
        
    # ---------------------------------------------------------- #
        
    N_RUNS = 16000

    perc_cells_to_update = 0.25


    # LIST OF TOLERANCES
    tol_f = 0.01 # tol% hange in Fano means convergence
    
    
    
    
    for time in range(N_RUNS):

        if time != 0:
            tmax_dict = {}
            for cell in range(N_CELLS):
                timepoints = sorted(cell_track[cell].trajectory.keys())
                tmax = max(timepoints)
                tmax_dict[cell] = tmax
            slow_cells= [k for v,k in sorted([(v,k) for k,v in tmax_dict.items()]) ][: int(round(perc_cells_to_update*N_CELLS))]

        else: 
            slow_cells = range(N_CELLS)


        for index in slow_cells:

            if cell_status[index] >= cell_track[index].deg_size:
                Propensity_deg = cell_track[index].k_deg*(math.factorial(cell_status[index])/(math.factorial(cell_status[index] - cell_track[index].deg_size)))

            else:
                Propensity_deg = 0.0


            # NORMALISE the propensities such that they sum to 1

            P_deg = Propensity_deg/(cell_track[index].k_burst + Propensity_deg)

            P_burst = cell_track[index].k_burst/(cell_track[index].k_burst + Propensity_deg)

            # determine step-size in time (first RANDOM number)

            lambdA = (cell_track[index].k_burst + Propensity_deg)
            delta_time = random.expovariate(lambdA)
            current_time = last_time_point[index]
            new_time = current_time + delta_time
            last_time_point[index] = new_time

            # pick a second RANDOM number that decides which event takes place: bursting or degradation

            r = random.random()

            if r <= P_burst: 
                cell_status[index] = cell_status[index] + cell_track[index].burst_size

            elif r > 1.0 - P_deg:
                cell_status[index] = cell_status[index] - cell_track[index].deg_size


            cell_track[index].trajectory[new_time] = cell_status[index] ### HERE is where I update the cell mRNA content!!


        # EACH RUN find the cell with the minimal simulated time

        list_of_max_times = []

        for i in range(N_CELLS):
            list_of_max_times.append(max(sorted(cell_track[i].trajectory.keys())))

        min_time = min(list_of_max_times)


        timepoints_on_left = []
        distribution_left = []
        index_timepoints_on_left = []   



        if time != 0 and time%10 == 0:

            for cell in range(N_CELLS):

                timepoints = sorted(cell_track[cell].trajectory.keys())
                timepoints_left = [i for i in timepoints if i <= min_time ] # all timepoints left of min_time

                this_timepoint_on_left = np.max(timepoints_left) # the biggest timepoint on the left
                index_timepoints_on_left.append(len(timepoints_left)) # count the number of points on the left
                timepoints_on_left.append(this_timepoint_on_left) 
                distribution_left.append(cell_track[cell].trajectory[timepoints_on_left[cell]])


                [cell_track[cell].trajectory.pop(x, None) for x in timepoints[:timepoints.index(this_timepoint_on_left)] ]


            v= np.var(distribution_left)
            m= np.mean(distribution_left)      
            f= v/m
            s= st.skew(distribution_left)
            k= st.kurtosis(distribution_left)

            Dict_min_time_Fano_var_mean[min_time] = f,v,m,s,k

            # the list "distribution_left" contains the mRNA content of each cell, 

            # this two list now are in the same dictionary
            sorted_Dict_min_time_Fano_var_mean = sorted(Dict_min_time_Fano_var_mean.items(), key= operator.itemgetter(0))

            fano_factor_list = [x[1][0] for x in sorted_Dict_min_time_Fano_var_mean]
            skewness_list = [x[1][3] for x in sorted_Dict_min_time_Fano_var_mean]



            if len(fano_factor_list) > 5:
                perc_fano = abs(fano_factor_list[-1] - fano_factor_list[-2])/fano_factor_list[-2]
                perc_fano2 = abs(fano_factor_list[-1] - fano_factor_list[-3])/fano_factor_list[-3]
                perc_fano3 = abs(fano_factor_list[-1] - fano_factor_list[-4])/fano_factor_list[-4]
                perc_fano4 = abs(fano_factor_list[-1] - fano_factor_list[-5])/fano_factor_list[-5]



                current_moments = [cell_track[0].k_deg,cell_track[0].k_burst,cell_track[0].deg_size,cell_track[0].burst_size,cell_track[0].Half_life, min_time,cell_track[0].macroscopic_mean, round(m,2),round(v,2),cell_track[0].van_Kampen_Fano_factor, round(f,2)]
                current_table = pd.DataFrame([current_moments], columns=['k_deg','k_burst','deg_size','burst_size','Theor. Half_life','Time','Theor. MEAN','Mean', 'Variance', 'Van Kampen FF', 'Fano Factor'] )   
                
                if all([x <= tol_f for x in [perc_fano,perc_fano2,perc_fano3,perc_fano4]]):

                    ss_distr = [int(x) for x in  distribution_left] ## if you want to plot this values need to be transformed in integer

                    break


    return ss_distr

