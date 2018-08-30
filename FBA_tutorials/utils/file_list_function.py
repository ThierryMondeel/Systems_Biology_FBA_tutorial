
# coding: utf-8

# In[ ]:

def file_list(directory = '.', extension ='_FPKM.json' ):
    
    import os
    
    f_list = [directory+f for f in os.listdir(directory) if extension in f ]
    
    return f_list   

