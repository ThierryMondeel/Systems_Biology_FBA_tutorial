
# coding: utf-8

# In[6]:

def glossary(key_word):
    """
    Generate a string that is the definition of the key word based on the fallowing **references**:
    
    Kærn et al._2005_Stochasticity in gene expression From theories to phenotypes.pdf
    Lenstra et al._2016_Transcription Dynamics in Living Cells
    Raj, Oudenaarden_2008_Review Nature , Nurture , or Chance Stochastic Gene Expression and Its Consequences
    Symmons, Raj_2016_Review What ’ s Luck Got to Do with It Single Cells , Multiple Fates , and Biological Nondeterminism



    Arguments
    ---------
    key_word : string of words written in bold through this tutorial

    Returns
    -------
    A string, definition of the the key word.
    It return the value of the of the key word in the dictionary called glossary_dict
    """
    
    glossary_dict = {}

    glossary_dict['isogenic']                  = 'Genetically identical.     Individual cells within an isogenic population are typically the progeny of a single ancestor.'
    glossary_dict['nucleosome']                = 'The fundamental unit into which DNA and histones are packaged in eukaryotic cells.    It is the basic structural subunit of chromatin and consists of 200 bp of DNA and an octamer of histone proteins.'
    glossary_dict['Fano factor']               = 'Mathematically defined as the variance of a distribution divided by the mean.    It is named after Ugo Fano, an Italian American physicist.'
    glossary_dict['steady-state distribution'] = 'the distribution of mRNA per cell across a population that     is equilibrated in the sense that the distribution will not change over time'
    glossary_dict['MS2']                       = 'a bacteriophage whose coat protein binds strongly with a particular RNA hairpin'
	glossary_dict['Brownian motion']           = 'Brownian motion or pedesis (from Ancient Greek: πήδησις /pέːdεːsis/ "leaping") is the random motion of particles suspended in a fluid (a liquid or a gas) resulting from their collision with the fast-moving molecules in the fluid.'

    return glossary_dict[key_word]
    

