
# coding: utf-8

# In[ ]:

def flux_pie_plot(df):
    
    import pandas as pd

    import numpy as np

    from math import log, sqrt
    from collections import OrderedDict

    from bokeh.plotting import figure, show, output_notebook
    from bokeh.palettes import brewer
    import warnings

    # Change name of subSystem in something shorter
    
    for i in range(len(df.subSystems.values)):
        if 'Glycolysis/gluconeogenesis' in df.subSystems.values[i]:
            df.at[i, 'subSystems'] = 'Glycolysis'
        if 'Squalene and cholesterol synthesis' in df.subSystems.values[i]:
            df.at[i, 'subSystems'] = 'Cholesterol synthesis'

    cell_line_color = OrderedDict([
        ("MCF7",   'black'),
        ("MCF7_T", 'crimson'),
        ("MCF7_F", 'gold'),
         ("LTED",   'blueviolet'),
    ])

    # In this specification I assume no more than 2 subSystem types
    ss_type_color = {}
    ss_c_list = ['lightskyblue', 'lightsalmon', ]

    for i in range(len(df.ss_type.unique())):
        ss_type_color[df.ss_type.unique()[i]] = ss_c_list[i] 
    ss_type_color




    width = 800
    height = 800
    inner_radius = 90
    outer_radius = 300 - 10

    minr = sqrt(log(1000 * 1E4))
    maxr = sqrt(log(0.001 * 1E4))
    a = (outer_radius - inner_radius) / (minr - maxr)
    b = inner_radius - a * maxr

    def rad(mic):
        return a * np.sqrt(np.log(mic * 1E4)) + b

    big_angle = 2.0 * np.pi / (len(df) + 1)
    small_angle = big_angle / 7

    #_______________________________________________________________________________

    p = figure(plot_width=width, plot_height=height, title="",
        x_axis_type=None, y_axis_type=None,
        x_range=(-500, 500), y_range=(-500, 500),
        min_border=0, outline_line_color="white",
        background_fill_color="white")

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None

    #_______________________________________________________________________________

    # # subSystem TYPE ('Central Carbon Metabolism' or 'Peripheral metabolism')

    angles = np.pi/2 - big_angle/2 - df.index.to_series()*big_angle
    colors = [ss_type_color[ss_type] for ss_type in df.ss_type]
    p.annular_wedge(
        0, 0, inner_radius, outer_radius, -big_angle+angles, angles, color=colors,
    )

    #_______________________________________________________________________________
    # subsystems BARS small wedges (bar plots) 

    p.annular_wedge(0, 0, inner_radius, rad(df.LTED),
                    -big_angle+angles+5*small_angle, -big_angle+angles+6*small_angle,
                    color=cell_line_color['MCF7_F'])
    p.annular_wedge(0, 0, inner_radius, rad(df.MCF7),
                    -big_angle+angles+3.5*small_angle, -big_angle+angles+4.5*small_angle,
                    color=cell_line_color['MCF7'])
    p.annular_wedge(0, 0, inner_radius, rad(df.MCF7_T),
                    -big_angle+angles+2*small_angle, -big_angle+angles+3*small_angle,
                    color=cell_line_color['MCF7_T'])
    p.annular_wedge(0, 0, inner_radius, rad(df.LTED),
                    -big_angle+angles+0.5*small_angle, -big_angle+angles+1.5*small_angle,
                    color=cell_line_color['LTED'])

    #_______________________________________________________________________________


    # circular axes 
    labels = np.power(10.0, np.arange(-3, 4))
    radii = a * np.sqrt(np.log(labels * 1E4)) + b

    p.circle(0, 0, radius=radii, fill_color=None, line_color="white")
    #  y-axis labels
    p.text(0, radii[:-1], [str(r) for r in labels[:-1]],
           text_font_size="10pt", text_align="center", text_baseline="middle")

    # radial axes
    p.annular_wedge(0, 0, inner_radius-10, outer_radius+10,
                    -big_angle+angles, -big_angle+angles, color="black")


    #_______________________________________________________________________________

    minr_i = sqrt(log(0.001 * 1E4))
    maxr_i = sqrt(log(1000 * 1E4))
    a_i = (outer_radius - inner_radius) / (minr_i - maxr_i)
    b_i = inner_radius - a_i * maxr_i



    big_angle_i = 2.0 * np.pi / (len(df) + 1)
    small_angle_i = big_angle / 7

    radii = a_i * np.sqrt(np.log(labels * 1E4)) + b_i


    # subSystem labels
    xr = radii[0]*np.cos(np.array(-big_angle/2 + angles))
    yr = radii[0]*np.sin(np.array(-big_angle/2 + angles))
    label_angle=np.array(-big_angle/2+angles)
    label_angle[label_angle < -np.pi/2] += np.pi # easier to read labels on the left side

    for i in range(len(xr)):
        if xr[i] > 0:
            p.text(xr[i], yr[i], pd.Series(df.subSystems.loc[i]), angle=label_angle[i],
           text_font_size="10pt", text_align="left", text_baseline="middle")
        else:
            p.text(xr[i], yr[i], pd.Series(df.subSystems.loc[i]), angle=label_angle[i],
           text_font_size="10pt", text_align="right", text_baseline="middle")

    #_______________________________________________________________________________
    # LEGEND

    # subSystem type 
    p.circle([+140, +140], [+350, +370], color=list(ss_type_color.values()), radius=5)
    p.text([+160, +160], [+350, +370], text=[gr for gr in ss_type_color.keys()],
           text_font_size="10pt", text_align="left", text_baseline="middle")

    # cell lines
    p.rect([-40, -40, -40,-40], [30, 10, -10,-30], width=30, height=13,
           color=list(cell_line_color.values()))
    p.text([-15, -15, -15, -15], [30, 10, -10,-30], text=list(cell_line_color),
           text_font_size="9pt", text_align="left", text_baseline="middle")

    #_______________________________________________________________________________

    warnings.filterwarnings("ignore")
    output_notebook()

    return show(p)

