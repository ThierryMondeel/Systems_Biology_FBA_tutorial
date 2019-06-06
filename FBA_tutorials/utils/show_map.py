import escher,escher.urls,json,os
from IPython.display import HTML


def show_map(sol,map_loc,color=0):
    ''' Returns an escher Builder object for solution 'sol', map 'map_loc' and the supplied color scheme.
        sol:        the solution object containing the simulation results.
        map_loc:    filename of the map json
        color:      color scheme to use
    '''

    if color == 0:
        colors = [{'type': 'min', 'color': '#cccccc', 'size': 5},# grey to green to orange
                  {'type': 'mean', 'color': '#007F00', 'size': 10},
                  {'type': 'max', 'color': '#f0a900', 'size': 15}]
    else:
        print('Color scheme not defined!')
        return

    if type(sol) != dict:
        try:
            d = sol.fluxes
        except:
            print('An empty solution was passed.')
            d = {}
    else:
        d = sol # shorthand

    d2 = {} # for some reason my types (from float 64 to float) were not updating and with a new dictionary they do
    for key in d.keys(): # remove output like this: 1.653e-15
        d2[key] = round(float(d[key]),6)


    network = escher.Builder(map_json=map_loc, reaction_data=d2, 
                       reaction_styles=['color', 'size', 'abs', 'text'],
                       # change the default colors, blue to purple to red
                       reaction_scale=colors,
                       hide_secondary_metabolites=False,secondary_metabolite_radius=10,
                       highlight_missing=True) 
    return network