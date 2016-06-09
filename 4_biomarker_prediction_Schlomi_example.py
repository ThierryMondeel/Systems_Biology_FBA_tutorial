
# coding: utf-8

# # Reproducing the example in Figure 1 of Shlomi et al. 2009
# This code reproduces the network and biomarker predictions in the example from (Schlomi,2009).
# You can read the paper here: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2683725/ 
# 
# To follow this tutorial be sure you have read the following parts of paper:
# - abstract
# - We reproduced essential parts of the method description in this notebook below. 
# 
# We will use cobrapy itself for this without the overarching framework of cameo. Many of the commands will be similar or the same since Cameo is built on top of Cobrapy.

# In[1]:

import cobra
from cobra import Model, Reaction, Metabolite
from cobrapyTools import * # collection of my own scripts
import pandas as pd


# We already made the example network shown in Figure A from the original publication below. 

# In[2]:

model = cobra.io.read_sbml_model("models/Shlomi_example.xml")


# ### The method
# ** This is a selective copy/paste from the original paper: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2683725/**
# 
# We present a new computational approach for systematically predicting the pattern of metabolic biomarkers characterizing each metabolic disorder whose causative gene is included in the human metabolic network model (Duarte et al, 2007). Let a boundary metabolite denote a metabolite that is known to be taken-up or secreted between the intracellular and extracellular compartments (as indicated in the network model). Let an exchange interval denote a possible range of uptake and secretion fluxes of a given boundary exchange interval. For each metabolic disorder and each boundary metabolite, we predict its exchange interval between human tissues and biofluids, for both healthy and disease cases (Materials and methods). This exchange interval is computed through a CBM method called flux variability analysis (FVA) (Mahadevan and Schilling, 2003), which accounts for the entire space of feasible flux states that satisfy mass-balance stoichiometric constraints and reaction directionality constraints (embedded in the model of (Duarte et al, 2007)). For the healthy case, the exchange interval is computed while the reactions affected by the disease are constrained to be active, whereas for the disease case, they are constrained to be inactive. By comparing the predicted exchange interval between the healthy state and the disease state for each boundary metabolite, one can determine whether the pertaining boundary metabolite concentration in biofluids (termed biomarker) is expected to be elevated, reduced or unchanged (see Materials and methods). If the predicted changes are marked such that there is no overlap between the exchange intervals of the healthy case and the disease case, the predicted biomarker change is considered to be highly confident.
# 
# An illustrative example of the predicted biomarker changes' ranges and their underlying rationale for the healthy state and some disease state is depicted in Figures 1A and B. The predicted exchange intervals of metabolite M1 (M2) suggest that its extracellular concentration is elevated (reduced) in the disease case. The disjoint exchange intervals obtained for the healthy case and the disease case for both M1 and M2 render these predictions as highly confident. The exchange intervals of metabolite M6 (M4) suggest that their extracellular concentrations are elevated (reduced) in the disease case. Examining, for example, the exchange interval of metabolite M6 shows that in the healthy case, M6 can be either taken-up from biofluids or secreted in a lower rate (as some of it is required in the healthy state; Supplementary Figure 1). In the disease case, M6 (synthesized through M5) can only be secreted to biofluids. It should be noted that mass-balance stoichiometric constraints that play an important role in determining the exchange intervals of different metabolites and are accounted for by the CBM method (and as will be shown, play an important role in determining biomarker changes in addition to the network topology) are not depicted in this kind of illustration.

# <img src="images/shlomi_example.png">

# ### Assignment
# Make sure you understand why, if flux variablity analysis predicts an exchange interval of metabolite X to become more positive,
# this "predicts" that the metabolite builds up in the biofluids! And vice versa.

# ### Assignment
# Make sure you agree that we have reconstructed the network in Figure 1A correctly. Write a for loop that prints all reactions in the model and check that they are correct.

# In[ ]:




# ### Reproducing figure B
# We have written a python function that implements the entire method proposed by Shlomi et al in the "findBiomarkers()" function. For now we will simply check to see if it produces the same qualitative results as the original paper. 

# We have reconstructed the network in Figure 1A. We will now block reaction 1 and show that we predict what is in Figure 1B. 

# In[3]:

exchanges = [ rxn for rxn in model.reactions if rxn.products == [] 
             and rxn.boundary == 'system_boundary']
exchangesIds = [rxn.id for rxn in exchanges]
T = findBiomarkers(model,[model.reactions.R1.id],exchangesIds,eps=1); T


# Note the difference between the data we produced and figure B regarding M7. We predict reduced extracellular levels while Shlomi et al predict no change. This is due to the fact the Shlomi et al put in a threshold of at least 10% change before they consider something a biomarker. We just show the raw data here.

# In[ ]:



