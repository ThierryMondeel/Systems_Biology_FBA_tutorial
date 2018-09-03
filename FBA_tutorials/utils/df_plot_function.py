
# coding: utf-8

# In[ ]:

def df_plot(target_ss, f_pattern):
    
    import pandas as pd
    
    df  = pd.DataFrame(target_ss,columns=['subSystems', 'ss_type'])
    for cell_line in f_pattern.keys():
        values = []
        for i in range(len(target_ss)):
            SUM = f_pattern[cell_line].fluxes[f_pattern[cell_line].subSystem == target_ss[i][0]].sum()

            values.append(SUM)

        df[cell_line] = values

    return df

