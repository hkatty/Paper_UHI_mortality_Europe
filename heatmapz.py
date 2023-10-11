# -*- coding: utf-8 -*-
"""
BSD 3-Clause License

Copyright (c) 2019, Drazen Zaric
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---------

From https://pypi.org/project/heatmapz
Details: https://towardsdatascience.com/better-heatmaps-and-correlation-matrix-plots-in-python-41445d0f2bec
Examples: https://colab.research.google.com/drive/1YSvER-U3cwGplSHyXwaCXYhOWfWO53Iy

With aesthetic modifications by Katty Huang noted by KH in code below

---------

heatmap(x, y, **kwargs)
Parameters:
x : A list, np.array or pandas.Series containing the values for the horizontal dimension
y : A list, np.array or pandas.Series containing the values for the vertical dimension

Optional parameters:
color : A list, np.array or pandas.Series containing values based on which to apply the heatmap color. Should have the same length as x and y.
palette : A list of colors to use as the heatmap palette. The values from color are mapped onto the palette so that min(color) -> palette[0] and max(color) -> palette[len(palette)-1], and the values in between are linearly interpolated. A good way to choose or create a palette is to simply use Seaborn palettes (https://seaborn.pydata.org/tutorial/color_palettes.html).
color_range : A tuple (color_min, color_max) that enables capping the values of color being mapped to palette. All color values less than color_min are capped to color_min, and all color values larger than color_max are capped to color_max. Then those values are mapped to palette as described under color.
size : A list, np.array or pandas.Series containing values based on which to apply the size to the shapes in the plot. Should have the same length as x and y.
size_range : A tuple (size_min, size_max) that enables capping the values of size being applied to the shapes in the plot. Essentially controls min and max size of the shapes.
size_scale : Used to scale the size of the shapes in the plot to make them fit the size of the fields in the matrix. Default value is 500. You will likely need to fiddle with this parameter in order to find the right value for your figure size and the size range applied.
x_order : Should contain all distinct values of x ordered in the way you want them shown on the x-axis from left to right.
y_order : Should contain all distinct values of y ordered in the way you want them shown on the y-axis from bottom to top.
marker : Specify the shape to use in the plot. It can be any of the matplotlib marker shapes (https://matplotlib.org/api/markers_api.html). The default is 's' for square.
xlabel : Label for the x-axis. Default is none.
ylabel : Label for the y-axis. Default is none.

corrplot(data, size_scale=500, marker='s')
A shorthand function for making correlation plots from pandas dataframes.
Parameters:
data : You should pass the result of calling df.corr() on a dataframe.
size_scale : Used to scale the size of the shapes in the plot to make them fit the size of the fields in the matrix. Default value is 500. You will likely need to fiddle with this parameter in order to find the right value for your figure size and the size range applied.
marker : Specify the shape to use in the plot. It can be any of the matplotlib marker shapes (https://matplotlib.org/api/markers_api.html). The default is 's' for square.

"""

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


def heatmap(x, y, **kwargs):
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = [1]*len(x)

    if 'palette' in kwargs:
        palette = kwargs['palette']
        n_colors = len(palette)
    else:
        n_colors = 256 # Use 256 colors for the diverging color palette
        palette = sns.color_palette("Blues", n_colors) 

    if 'color_range' in kwargs:
        color_min, color_max = kwargs['color_range']
    else:
        color_min, color_max = min(color), max(color) # Range of values that will be mapped to the palette, i.e. min and max possible correlation

    def value_to_color(val):
        if color_min == color_max:
            return palette[-1]
        else:
            val_position = float((val - color_min)) / (color_max - color_min) # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            ind = int(val_position * (n_colors - 1)) # target index in the color palette
            return palette[ind]

    if 'size' in kwargs:
        size = kwargs['size']
    else:
        size = [1]*len(x)

    if 'size_range' in kwargs:
        size_min, size_max = kwargs['size_range'][0], kwargs['size_range'][1]
    else:
        size_min, size_max = min(size), max(size)

    size_scale = kwargs.get('size_scale', 500)

    def value_to_size(val):
        if size_min == size_max:
            return 1 * size_scale
        else:
            val_position = (val - size_min) * 0.99 / (size_max - size_min) + 0.01 # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            return val_position * size_scale
    if 'x_order' in kwargs: 
        x_names = [t for t in kwargs['x_order']]
    else:
        x_names = [t for t in sorted(set([v for v in x]))]
    x_to_num = {p[1]:p[0] for p in enumerate(x_names)}

    if 'y_order' in kwargs: 
        y_names = [t for t in kwargs['y_order']]
    else:
        y_names = [t for t in sorted(set([v for v in y]))]
    y_to_num = {p[1]:p[0] for p in enumerate(y_names)}

    #KH>>create fig instance if not passed from call so can apply tight_layout at the end
    if 'fig' in kwargs:
        fig = kwargs['fig']
    else:
        fig = plt.figure()
    #<<KH
    plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1) # Setup a 1x10 grid
    ax = plt.subplot(plot_grid[:,:-1]) # Use the left 14/15ths of the grid for the main plot

    marker = kwargs.get('marker', 's')

    kwargs_pass_on = {k:v for k,v in kwargs.items() if k not in [
         'color', 'palette', 'color_range', 'size', 'size_range', 'size_scale', 'marker', 'x_order', 'y_order', 'xlabel', 'ylabel', 'fig', 'clabel', 'title' #KH added fig, clabel, and title
    ]}

    ax.scatter(
        x=[x_to_num[v] for v in x],
        y=[y_to_num[v] for v in y],
        marker=marker,
        s=[value_to_size(v) for v in size], 
        c=[value_to_color(v) for v in color],
        **kwargs_pass_on
    )
    ax.set_xticks([v for k,v in x_to_num.items()])
    ax.set_xticklabels([k for k in x_to_num], rotation=45, horizontalalignment='right')
    ax.set_yticks([v for k,v in y_to_num.items()])
    ax.set_yticklabels([k for k in y_to_num])

    ax.grid(False, 'major')
    # ax.grid(True, 'minor') #KH remove background grid
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    # ax.set_facecolor('#F1F1F1') #KH remove background colour

    ax.set_xlabel(kwargs.get('xlabel', ''))
    ax.set_ylabel(kwargs.get('ylabel', ''))
    ax.set_title(kwargs.get('title','')) #KH add optional title

    # Add color legend on the right side of the plot
    if color_min < color_max:
        ax = plt.subplot(plot_grid[:,-1]) # Use the rightmost column of the plot

        col_x = [0]*len(palette) # Fixed x coordinate for the bars
        bar_y=np.linspace(color_min, color_max, n_colors) # y coordinates for each of the n_colors bars

        bar_height = bar_y[1] - bar_y[0]
        ax.barh(
            y=bar_y,
            width=[5]*len(palette), # Make bars 5 units wide
            left=col_x, # Make bars start at 0
            height=bar_height,
            color=palette,
            linewidth=0
        )
        ax.set_xlim(1, 2) # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
        ax.set_ylim(bar_y[0],bar_y[-1]) #KH added to remove white space at ends of colourbar
        ax.grid(False) # Hide grid
        ax.set_facecolor('white') # Make background white
        ax.set_xticks([]) # Remove horizontal ticks
        ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3)) # Show vertical ticks for min, middle and max
        ax.yaxis.tick_right() # Show vertical ticks on the right 
        ax.set_ylabel(kwargs.get('clabel', ''),rotation=270) #KH add colourbar label is specified
        ax.yaxis.set_label_position("right") #KH label on the right
    
    plot_grid.tight_layout(fig) #KH


def corrplot(data, size_scale=500, marker='s'):
    corr = pd.melt(data.reset_index(), id_vars='index').replace(np.nan, 0)
    corr.columns = ['x', 'y', 'value']
    heatmap(
        corr['x'], corr['y'],
        color=corr['value'], color_range=[-1, 1],
        palette=sns.diverging_palette(20, 220, n=256),
        size=corr['value'].abs(), size_range=[0,1],
        marker=marker,
        x_order=data.columns,
        y_order=data.columns[::-1],
        size_scale=size_scale
    )

