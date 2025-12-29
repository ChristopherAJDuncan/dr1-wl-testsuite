# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2025-06-24 08:54:57
# @Last Modified by:   lshuns
# @Last Modified time: 2025-12-28 18:35:02

### With the given split criteria, compare key properties of subsamples
#### north vs south fields

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

import plotting

## >>>>>>>>>>>>> I/O

## which shear measurement
shear_type = 'lensmc'
# shear_type = 'metacal'

## number of tomo bins
Ntomo = 6

## where to find catalogues
TITLE = 'euclid_tr1_a_v1 + euclid_tr1_b_v1_1'
inpath_list = ['/sdf/data/kipac/u/liss/euclid/TR1/23538.parquet', 
               '/sdf/data/kipac/u/liss/euclid/TR1/23539.parquet']
used_cols_list = [['object_id', 
                   'right_ascension', 
                   'declination', 
                   'flux_radius',
                   'phz_pp_mode_mu', 
                   'phz_pp_mode_mr', 
                   'sersic_sersic_vis_radius', 
                   'sersic_sersic_vis_index', 
                   'sersic_sersic_vis_axis_ratio',
                   'flux_detection_total'],
                  ['object_id', 
                   'tom_bin_id', 
                   'phz_mode_1', 
                   f'she_{shear_type}_weight']]

## which columns to plot
plot_cols = ['flux_radius', 
             'phz_mode_1', 
             'sersic_sersic_vis_radius', 
             'sersic_sersic_vis_index', 
             'sersic_sersic_vis_axis_ratio',
            'mag_detection_total']
XRANGE_list = [[1, 20], 
               [0, 2.5],  
               [0.01, 20],
               [0, 9],
               [0, 1],
               [20, 25]]
xlog_list = [True, 
             False,
             True,
             False,
             False,
             False
             ]

## >>>>>>>>>>>>>> Workhorse

## +++ Load and merge catalogue
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
cata = cata[(cata['flux_detection_total']>0)].reset_index(drop=True)
print(f">>> Number after flux quality filter", len(cata))
## calculate magnitude
cata['mag_detection_total'] = 23.9 - 2.5 * np.log10(cata['flux_detection_total'].values)

## +++ Assign subsample ID
cata['subsample_ID'] = 0
## calculate colour and cut
cata.loc[(cata['declination'].values < 0), 'subsample_ID'] = 1
print(">>>>>>> subsample number ratio (north/south)", np.sum(cata['subsample_ID']==0)/np.sum(cata['subsample_ID']==1))
print(">>>>>>> subsample weight ratio (north/south)", np.sum(cata.loc[cata['subsample_ID']==0, f'she_{shear_type}_weight'].values)/np.sum(cata.loc[cata['subsample_ID']==1, f'she_{shear_type}_weight'].values))

## +++ 1D hist without tomo bins
mask_south = (cata['subsample_ID']==1)
for i_col, plot_col in enumerate(plot_cols):
    outpath = f'./plots/B_north_south_1Dhist_all_{plot_col}_{shear_type}_weighted.png'
    XLABEL = plot_col
    xlog = xlog_list[i_col]
    XRANGE = XRANGE_list[i_col]
    paras = [cata.loc[mask_south, plot_col].values, 
                cata.loc[(~mask_south), plot_col].values]
    wgs = [cata.loc[mask_south, f'she_{shear_type}_weight'].values, 
                cata.loc[(~mask_south), f'she_{shear_type}_weight'].values]
    COLORs = ['cyan', 'orange']
    LABELs = ['south', 'north']
    nbins = 60
    DENSITY = False
    YLABEL = f'{shear_type}-weighted counts'
    plotting.HistPlotFunc(outpath,
                paras, wgs, COLORs, LABELs,
                nbins, XRANGE, YRANGE=None,
                XLABEL=XLABEL, YLABEL=YLABEL,
                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                TITLE=TITLE, xtick_min_label=True, ytick_min_label=True,
                xtick_spe=None, ytick_spe=None,
                vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                xlog=xlog, ylog=False,
                loc_legend='best', 
                LABEL_position='inSub', LABEL_cols=1,
                font_size=12, usetex=False,
                cumulative=False, 
                FIGSIZE=[6.4, 4.8],
                LINEs=None, LINEWs=None,
                TIGHT=False,
                alpha=None,
                HISTTYPE_list=None)

## +++ 1D hist with tomo bins
N_plots = Ntomo
mask_south = (cata['subsample_ID']==1)
for i_col, plot_col in enumerate(plot_cols):
    outpath = f'./plots/B_north_south_1Dhist_tomo_{plot_col}_{shear_type}_weighted.png'
    XLABEL = plot_col
    xlog = xlog_list[i_col]
    XRANGE = XRANGE_list[i_col]
    paras_list = []
    wgs_list = []
    subLABEL_list = []
    for i_ZBbin in range(Ntomo):
        mask_tomo = cata['tom_bin_id']==(i_ZBbin+1)
        paras_list.append([cata.loc[mask_south&mask_tomo, plot_col].values, 
                           cata.loc[(~mask_south)&mask_tomo, plot_col].values])
        wgs_list.append([cata.loc[mask_south&mask_tomo, f'she_{shear_type}_weight'].values, 
                           cata.loc[(~mask_south)&mask_tomo, f'she_{shear_type}_weight'].values])
        subLABEL_list.append(f'tomo {i_ZBbin+1}')
    COLORs_list = [['cyan', 'orange']] * N_plots
    LABELs_list = [['south', 'north']] * N_plots
    nbins_list = [60] * N_plots
    DENSITY = False
    YLABEL = f'{shear_type}-weighted counts'
    plotting.HistPlotFunc_subplots(outpath, N_plots,
                                paras_list, wgs_list, COLORs_list, LABELs_list,
                                nbins_list, XRANGE, YRANGE=None,
                                LINEs_list=None, LINEWs_list=None,
                                subLABEL_list=subLABEL_list, subLABEL_locX=0.1, subLABEL_locY=0.8,
                                XLABEL=XLABEL, YLABEL=YLABEL,
                                DENSITY=DENSITY, HISTTYPE='step', STACKED=False,
                                TITLE=None, xtick_min_label=True, ytick_min_label=True,
                                xtick_spe=None, ytick_spe=None,
                                vlines_list=None, 
                                vline_styles_list=None, vline_colors_list=None, vline_labels_list=None, vline_widths_list=None,
                                hlines_list=None, 
                                hline_styles_list=None, hline_colors_list=None, hline_labels_list=None, hline_widths_list=None,
                                xlog=xlog, ylog=False,
                                loc_legend='best',
                                font_size=12, usetex=False,
                                LABEL_position='top', LABEL_position_SUBid=0,
                                LABEL_cols=2, 
                                FIGSIZE=[6.4, 4.8],
                                TIGHT=False,
                                HISTTYPEs_list=None,
                                shareX=True,
                                shareY=True,
                                XRANGE_list=None,
                                YRANGE_list=None)