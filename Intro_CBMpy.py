
# coding: utf-8

# In[24]:

import vonda
import cbmpy as cbm
from cobrapyTools import *


# ## CBMpy INTRO

# In[5]:

model = cbm.CBRead.readSBML3FBC('../SYNPOL/clostridium_ljungdahlii/models/nagarajan_2013/c_ljungdahlii_nagarajan_2013_update.xml')


# In[10]:

model.reactions[1].getId()


# In[17]:

sol = cbm.CBSolver.analyzeModel(model,return_lp_obj = True)


# In[21]:

cbm.CBCPLEX.cplx_MinimizeNumActiveFluxes(model)


# In[37]:

cbm.CBModelTools.printSolution(model)


# ## VONDA INTRO
# Standard vizualization with pre-made data

# In[25]:

vmod = vonda.PVisualizer('Synechocystis.svg')
D_fluxes = vmod.importKeyValueData("iTM686_FBA.txt")
vmod.doMapReactions(D_fluxes)
vmod.doMapReactions(D_fluxes,valuesRange=['higher',0.01],
    filename_out = "Synechocystis_reactions_high")
vmod.output_dir


# The last simulation done above

# In[39]:

vmod = vonda.PVisualizer('Synechocystis.svg')
D_fluxes = model.getReactionValues()
vmod.doMapReactions(D_fluxes)
vmod.doMapReactions(D_fluxes,valuesRange=['higher',0.01],
    filename_out = "Synechocystis_reactions_high")
vmod.output_dir


# In[36]:

model.getId()


# In[ ]:



