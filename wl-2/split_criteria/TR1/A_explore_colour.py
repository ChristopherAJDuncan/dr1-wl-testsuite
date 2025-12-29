# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2025-06-24 08:54:57
# @Last Modified by:   lshuns
# @Last Modified time: 2025-12-28 12:21:37

### Explore the restframe colour info with the TR1 data
###### to guide the split criteria

import numpy as np
import pandas as pd 

import plotting

## >>>>>>>>>>>>> I/O

## which shear measurement
shear_type = 'lensmc'
# shear_type = 'metacal'

## where to find catalogues
TITLE = 'euclid_tr1_a_v1 + euclid_tr1_b_v1_1'
inpath_list = ['/sdf/data/kipac/u/liss/euclid/TR1/23538.parquet', 
               '/sdf/data/kipac/u/liss/euclid/TR1/23539.parquet']
used_cols_list = [['object_id', 
                   'flux_detection_total', 
                   'fluxerr_detection_total',
                   'phz_pp_mode_mu', 
                   'phz_pp_mode_mr', ],
                  ['object_id', 
                   'phz_mode_1', 
                   f'she_{shear_type}_weight']]

## >>>>>>>>>>>>>> workhorse

## +++ Load and merge catalogues
for i_cata, inpath in enumerate(inpath_list):
    if i_cata == 0:
        cata = pd.read_parquet(inpath)[used_cols_list[i_cata]]
        print(f">>> Number loaded from cata {i_cata}", len(cata))
    else:
        cata_tmp = pd.read_parquet(inpath)[used_cols_list[i_cata]]
        print(f">>> Number loaded from cata {i_cata}", len(cata_tmp))
        cata = cata.merge(cata_tmp, on='object_id')
        del cata_tmp
        print(f">>>>>> Number after merging", len(cata))
## weight selection
cata = cata.loc[(cata[f'she_{shear_type}_weight']>0)].reset_index(drop=True)
print(f">>> Number after she_{shear_type}_weight selection", len(cata))
## additional data quality filters
cata = cata[(cata['flux_detection_total']>0) & (cata['fluxerr_detection_total']>0)].reset_index(drop=True)
print(f">>> Number after flux quality filter", len(cata))
## calculate restframe u-r colour
restframe_u_r = cata['phz_pp_mode_mu'].values - cata['phz_pp_mode_mr'].values

## +++ 1-D colour distributions
# outpath = 'show'
outpath = f'./plots/A_colour_1Dhist_u-r_{shear_type}_weighted.png'
paras = [restframe_u_r]
COLORs = ['k']
LABELs = None
nbins = 60
XRANGE = [-0.3, 3.]
XLABEL = 'Restframe u-r colour'
# wgs = None
# YLABEL = 'Number counts'
wgs = [cata[f'she_{shear_type}_weight'].values]
YLABEL = f'{shear_type}-weighted counts'
plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=False, HISTTYPE='step', STACKED=False,
                TITLE=TITLE, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=False, ylog=False,
                loc_legend='best', 
                LABEL_position='inSub', LABEL_cols=1,
                font_size=12, usetex=False,
                cumulative=False, 
                FIGSIZE=[6.4, 4.8],
                LINEs=None, LINEWs=None,
                TIGHT=False,
                alpha=None,
                HISTTYPE_list=None)

## +++ 2D hist of colour VS key parameters
wg = cata[f'she_{shear_type}_weight'].values
CBAR_LABEL = f'{shear_type}-weighted counts'
nbins = 60
count_log = True
## x-axis info
x_val_list = [cata['phz_mode_1'].values, 
              cata['phz_pp_mode_mr'].values,
              23.9 - 2.5 * np.log10(cata['flux_detection_total'].values),
              cata['flux_detection_total'].values / cata['fluxerr_detection_total'].values]
XLABEL_list = ['phz_mode_1', 
               'phz_pp_mode_mr',
               'VIS apparent mag',
               'VIS flux SNR']
XRANGE_list = [[0.2, 2.5], 
               [-25, -15],
               [20, 25],
               [0, 100]]
## y-axis info
y_val_list = [restframe_u_r]
YLABEL_list = ['Restframe u-r colour']
YRANGE_list = [[-0.3, 3.0]]
## loop over and plot
for i_x, x_val in enumerate(x_val_list):
    XLABEL = XLABEL_list[i_x]
    XRANGE = XRANGE_list[i_x]
    for i_y, y_val in enumerate(y_val_list):
        YLABEL = YLABEL_list[i_y]
        YRANGE = YRANGE_list[i_y]
        outpath = f'./plots/A_colour_2Dhist_u-r_vs_{XLABEL.replace(" ", "_")}_{shear_type}_weighted.png'
        plotting.Hist2DPlotFunc(outpath,
                        x_val, y_val, wg,
                        nbins, XRANGE=XRANGE, YRANGE=YRANGE,
                        XLABEL=XLABEL, YLABEL=YLABEL, CBAR_LABEL=CBAR_LABEL,
                        COLOR_MAP='Reds',
                        DENSITY=False, count_scale=[None, None], count_log=count_log,
                        TITLE=TITLE, xtick_min_label=True, ytick_min_label=True,
                        xtick_spe=None, ytick_spe=None,
                        vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                        hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                        font_size=12, usetex=False, 
                        FIGSIZE=[6.4, 4.8],
                        TIGHT=False,
                        xlog=False, ylog=False)
    
## +++ Scatter plots coloured by lensmc weights
## randomly pick 1e5 points
random_ids = np.random.choice(range(len(cata)), 
                              100000, replace=False)
POINTS = 0.1
cval = cata[f'she_{shear_type}_weight'].values[random_ids]
bar_label = f'{shear_type} weight'
## x-axis info
x_val_list = [cata['phz_mode_1'].values[random_ids], 
              cata['phz_pp_mode_mr'].values[random_ids],
              23.9 - 2.5 * np.log10(cata['flux_detection_total'].values[random_ids]),
              cata['flux_detection_total'].values[random_ids] / cata['fluxerr_detection_total'].values[random_ids]]
XLABEL_list = ['phz_mode_1', 
               'phz_pp_mode_mr',
               'VIS apparent mag',
               'VIS flux SNR']
XRANGE_list = [[0.2, 2.5], 
               [-25, -15],
               [20, 25],
               [0, 100]]
## y-axis info
y_val_list = [restframe_u_r[random_ids]]
YLABEL_list = ['Restframe u-r colour']
YRANGE_list = [[-0.3, 3.0]]
## loop over and plot
for i_x, xval in enumerate(x_val_list):
    XRANGE = XRANGE_list[i_x]
    XLABEL = XLABEL_list[i_x]
    for i_y, yval in enumerate(y_val_list):
        YRANGE = YRANGE_list[i_y]
        YLABEL = YLABEL_list[i_y]
        outpath = f'./plots/A_colour_scatter_u-r_vs_{XLABEL.replace(" ", "_")}_coloured_by_{shear_type}_weight.png'
        plotting.ScatterPlotFunc(outpath,
                        xval, yval, POINT=None, POINTS=POINTS, alpha=None,
                        cval=cval, cmap=None, cmin=None, cmax=None, clog=False,
                        bar_loc=None, bar_ori=None, bar_label=bar_label, bar_tick=None,
                        XRANGE=XRANGE, YRANGE=YRANGE,
                        XLABEL=XLABEL, YLABEL=YLABEL, TITLE=TITLE,
                        xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                        vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                        hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                        xlog=False, invertX=False, ylog=False, invertY=False, 
                        loc_legend='best', legend_frame=False,
                        font_size=12, usetex=False,
                        FIGSIZE=[6.4, 4.8],
                        texPacks=None)