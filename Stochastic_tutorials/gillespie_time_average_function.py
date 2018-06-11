
# coding: utf-8

# In[52]:

def gillespie_time_average(N_RUNS = 5000, k_burst = 50.0, k_deg = 1.0, burst_size = 1.0):
    
    import random
    import numpy as np
    import math
    import operator
    
    # Parameters
    k_burst    = k_burst
    k_deg      = k_deg

    burst_size = burst_size
    deg_size     = 1

    # Initial conditions
    t_0 = 0
    initial_mrna_zero = 0
    
    cell_status = {}
    cell_status[t_0] = initial_mrna_zero

    N_RUNS = N_RUNS

    #1 : Initialise species and parameters!

    cell_status ={}
    t_0 = 0
    initial_mrna_zero = 0

    cell_status[t_0] = initial_mrna_zero

    #4 : Iterations of step 2 & 3

    for run in range(N_RUNS):

        last_time_point = max(cell_status.keys())


        if cell_status[last_time_point] >= deg_size:
            Propensity_deg = k_deg*(math.factorial(cell_status[last_time_point])/(math.factorial(cell_status[last_time_point] - deg_size)))

        else:
            Propensity_deg = 0.0

        Propensity_burst = k_burst

        SUM_propensities = Propensity_deg + Propensity_burst

        # NORMALISE the propensities such that they sum to 1

        P_deg = Propensity_deg/SUM_propensities

        P_burst = Propensity_burst/SUM_propensities

        # determine step-size in time (first RANDOM number)

        lambdA = (SUM_propensities)
        delta_time = random.expovariate(lambdA)

        current_time = last_time_point
        new_time = round((current_time + delta_time), 3)

        # pick a second RANDOM number that decides which event takes place: bursting or degradation

        r = random.random()

        # HERE is where I update the cell mRNA content in case of BURSTING

        if r <= P_burst: 
            cell_status[new_time] = cell_status[last_time_point] + burst_size 
        # HERE is where I update the cell mRNA content in case of DEGRADATION
        elif r > 1.0 - P_deg: 
            cell_status[new_time] = cell_status[last_time_point] - deg_size 
    
    sorted_cell_status = sorted(cell_status.items(), key= operator.itemgetter(0))
    
    time_steps = [x[0] for x in sorted_cell_status]
    mRNA_trajectory = [x[1] for x in sorted_cell_status]
    # the output of this function are 2 lists:
    return time_steps, mRNA_trajectory

