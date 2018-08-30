
# coding: utf-8

# In[7]:

def flux_pattern(file_list, analysis = 'FBA',fraction_of_optimum=1.0):
    ''' 
Parameters
----------
file list : list of genome-scale metabolic models in json format

analysis type : 

    - FBA  : Flux Balance Analisys, this will maximize or minimize (maximizing is the default) 
             flux through the objective reactions;
             
    - pFBA : Parsimonius Flux Balance Analisys finds a flux distribution which gives the optimal growth rate, 
             but minimizes the total sum of flux (for more details on pFBA, please see Lewis et al.,2010).
    
    - FVA  : Flux Variability Analisys finds the ranges of each metabolic flux at the optimum 
             (the default parameter is fraction_of_optimum=1.0, can be set to a different fraction).


OUTPUT 
----------
    resultDict = { cell_line: DataFrame of flux pattern}
   
   '''
    
    import pandas as pd
    import cobra
    import warnings
    from cobra.flux_analysis import flux_variability_analysis 
    from colorama import Fore
    
    # load models and save store them in a dictionary
    
    # to avoid warnings!
    warnings.filterwarnings("ignore")
    
    print()
    print('Loading the following models:'+ str(file_list).replace('[','').replace(']',''))
    print()
    
    # Make a Model's Dictionary
    model_dict = {}
    
    for model in file_list:

        model_name = model.replace('_FPKM.json','')
        model_name = model_name.split('/')[-1]

        model_dict[model_name] = cobra.io.load_json_model(model)
    
    print('All models are loaded')
    print()

# choose which Analysis to perform
   
    print('_______________________________________________________________________________ ')
    print()
    print('The flux pattern will be computed using '+ Fore.RED + analysis)
    print(Fore.BLACK + '_______________________________________________________________________________ ')

#-------------------------------------------------------------------------------- 

    if analysis == 'FBA':

     # Dictionary of flux pattern in each cell-line
        f_pattern = {}

        for m in model_dict.keys():

            # Optimize and save flux pattern

            df_fluxes = model_dict[m].optimize().fluxes.to_frame()

            # subDataframe of reactions carrying flux
            rxns_flux = df_fluxes[df_fluxes.fluxes != 0]

            # Add a new column for subsystem

            rxns_flux['subSystem'] = ['']*len(rxns_flux)

            # Update subsystem column, and fluxes with an absolute number

            for rxn in rxns_flux.index:
                ss = model_dict[m].reactions.get_by_id(rxn).subsystem
                rxns_flux.at[rxn,'subSystem'] = ss[0]
                rxns_flux.at[rxn,'fluxes'] = abs(rxns_flux.loc[rxn]['fluxes'])

            rxns_flux = rxns_flux.sort_values('fluxes', ascending=False)
            print()
            print('Flux pattern of '+ 'cell-line ' + m + ' is ready!')
            print()
            print(rxns_flux.head())
            print()
            print('--------------------------------------------------------------------------------')

            f_pattern[m] = rxns_flux

#-------------------------------------------------------------------------------- 

    elif analysis == 'pFBA':
        
        pfba = cobra.flux_analysis.pfba
        
     # Dictionary of flux pattern in each cell-line
        f_pattern = {}

        for m in model_dict.keys():

            # Optimize and save flux pattern

            df_fluxes = pfba(model_dict[m]).fluxes.to_frame()

            # subDataframe of reactions carrying flux
            rxns_flux = df_fluxes[df_fluxes.fluxes != 0]
            df_fluxes = model_dict[m].optimize().fluxes.to_frame()
            
            # Add a new column for subsystem

            rxns_flux['subSystem'] = ['']*len(rxns_flux)
            # Update subsystem column, and fluxes with an absolute number

            for rxn in rxns_flux.index:
                ss = model_dict[m].reactions.get_by_id(rxn).subsystem
                rxns_flux.at[rxn,'subSystem'] = ss[0]
                rxns_flux.at[rxn,'fluxes'] = abs(rxns_flux.loc[rxn]['fluxes'])

            rxns_flux = rxns_flux.sort_values('fluxes', ascending=False)
            print()
            print('Flux pattern of '+ 'cell-line ' + m + ' is ready!')
            print()
            print(rxns_flux.head())
            print()
            print('_______________________________________________________________________________ ')

            f_pattern[m] = rxns_flux



 #-------------------------------------------------------------------------------- 
    elif analysis == 'FVA':
                 
        fva = cobra.flux_analysis.flux_variability_analysis
    
    # Dictionary of flux pattern in each cell-line
        f_pattern = {}

        for m in model_dict.keys():

            # Perform FVA to compute min and flux values carried by each single rxn

            

            df_fluxes = fva(model_dict[m], 
                            model_dict[m].reactions, 
                            fraction_of_optimum = fraction_of_optimum)

            # I compute the max absolute value between min and max flux resulting from the FVA

            maxCol=lambda x: max(x.min(), x.max(), key=abs)
            df_fluxes['fluxes']  = df_fluxes.apply(maxCol,axis=1)
            df_fluxes = df_fluxes.iloc[:,2:] # the third column is "abs_max_flux"


            # subDataframe of reactions carrying flux

            rxns_flux = df_fluxes[df_fluxes.fluxes != 0]
            
            # Add a new column for subsystem

            rxns_flux['subSystem'] = ['']*len(rxns_flux)

            # Update subsystem column, and fluxes with an absolute number

            for rxn in rxns_flux.index:
                ss = model_dict[m].reactions.get_by_id(rxn).subsystem
                if len(rxns_flux.loc[rxn]['subSystem']) == 0:
                    rxns_flux.at[rxn,'subSystem'] = ''
                    rxns_flux.at[rxn,'fluxes'] = abs(rxns_flux.loc[rxn]['fluxes'])
            print()
            print('Flux pattern of '+ 'cell-line ' + m + ' is ready!')
            print()
            print(rxns_flux.head())
            print()
            print('_______________________________________________________________________________ ')

            f_pattern[m] = rxns_flux 
   
    
    
    
    
    
    return f_pattern

