
# coding: utf-8

# ### Applying Shlomi et al's biomarker prediction method to PKU on RECON2
# The original publication looked at RECON2's predecessor RECON1. We will reproduce their analysis on RECON2 instead.

# In[ ]:

import cobra
from cobra import Model, Reaction, Metabolite
from cobrapyTools import * # collection of my own scripts
import pandas as pd
pd.set_option('display.max_colwidth', -1)

model = cobra.io.read_sbml_model("models/Recon2.v04_pythonComp.xml")


# ### Assignment: biomarker prediction do it yourself
# Take the previous tutorial's last command cell as a template (we already copy-pasted it for you) to do the same analysis here on the full RECON2 model and the PKU disease state. 
# 
# **Tips**
# - Don't reinvent the wheel, the idea is the same as the last tutorial.
# - Instead of giving 'R1' as the disease reaction you should now give the PKU enzyme reactions. 
# - Also think about the fact that there are two equivalent reactions that you have to account for not just one.
# - The number of exchange reactions here is much bigger than in the example. The FVA computation will take quite a bit longer as a result. It may take 2 minutes or so. 

# In[ ]:

exchanges = [ rxn for rxn in model.reactions if rxn.products == [] 
             and rxn.boundary == 'system_boundary']
for rxn in exchanges:
    if rxn.lower_bound < 0:
        rxn.lower_bound = -999999
    if rxn.upper_bound > 0:
        rxn.upper_bound = 999999
        
exchangesIds = [rxn.id for rxn in exchanges]
T = findBiomarkers(model,[??Do something here??],exchangesIds,eps=1); T


# ### Question
# If you see the prediction of the biofluids/tissue as the brain tissue: does the model correctly predict issues with neurotransmitters in the brain? 

# ### Bonus assignment
# What biomarkers are predicted when you focus on blocking the cofactor biopterin recycling reactions that also produce PKU? To get you started we included some code below.

# In[ ]:

print model.reactions.DHPR.genes
print model.reactions.DHPR2.genes
print model.reactions.r0398.genes

T = findBiomarkers(model,[??Do something here??],exchangesIds,eps=1); T


# ### Question
# Are there any differences between the predictions for the two different ways to get PKU?

# In[ ]:



