# -*- coding: utf-8 -*-
# @Author: lshuns
# @Date:   2025-06-24 08:40:02
# @Last Modified by:   lshuns
# @Last Modified time: 2025-12-28 18:34:51

### Measure 2pt for the whole sample and defined subsamples
###### red vs blue galaxies

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Custom modules
sys.path.insert(0, '../utils')
from CorrFunc import CalCorrFunc, CorrPlotFunc_half, CorrPlotFunc_half_diff2base

# >>>>>>>>>>>> I/O

## which shear measurement
shear_type = 'lensmc'
# shear_type = 'metacal'

## where to find catalogues
TITLE = 'euclid_tr1_a_v1 + euclid_tr1_b_v1_1'
inpath_list = ['/sdf/data/kipac/u/liss/euclid/TR1/23538.parquet', 
               '/sdf/data/kipac/u/liss/euclid/TR1/23539.parquet']
used_cols_list = [['object_id', 
                   'right_ascension', 
                   'declination', 
                   'phz_pp_mode_mu', 
                   'phz_pp_mode_mr'],
                  ['object_id', 
                   'tom_bin_id', 
                   'phz_mode_1', 
                   f'she_{shear_type}_weight',
                   f'she_{shear_type}_e1_corrected',
                   f'she_{shear_type}_e2_corrected']]

## where to save
outdir = f'/sdf/data/kipac/u/liss/euclid/TR1/XI'
Path(outdir).mkdir(parents=True, exist_ok=True)

## u-r cut value (selection criteria)
cut_value = 1.2

## number of tomo bins
Ntomo = 6

## parameters for treecorr
theta_nbins = 20
theta_min = 1
theta_max = 400
theta_unit = "arcmin"
theta_bin_slop = 0.05
nthr = 128

## >>>>>>>>>>>> Workhorse

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
print(f">>> Effective number", (np.sum(cata[f'she_{shear_type}_weight']))**2 / np.sum(cata[f'she_{shear_type}_weight']**2))
## rename the columns
cols_rename = {f'she_{shear_type}_weight': 'shape_weight', 
                f'she_{shear_type}_e1_corrected': 'e1',
                f'she_{shear_type}_e2_corrected': 'e2',
                'right_ascension': 'RA',
                'declination': 'DEC'}
cata.rename(columns=cols_rename, inplace=True)

## +++ Assign subsample ID
cata['subsample_ID'] = 0
## calculate colour and cut
cata.loc[((cata['phz_pp_mode_mu'].values - cata['phz_pp_mode_mr'].values)>=cut_value), 'subsample_ID'] = 1
print(f">>> Subsamples defined with cut_value = {cut_value}")
print(">>>>>>> Subsample number ratio (blue/red)", np.sum(cata['subsample_ID']==0)/np.sum(cata['subsample_ID']==1))
print(">>>>>>> Subsample weight ratio (blue/red)", np.sum(cata.loc[cata['subsample_ID']==0, 'shape_weight'].values)/np.sum(cata.loc[cata['subsample_ID']==1, 'shape_weight'].values))

## +++ Calculate correlation function
print('+++ Total number of tomo bins', Ntomo)
outpath_2pt_list = []
## whole sample
print('++++++ Measuring 2pt for the whole sample')
cata_final = []
for i_ZBbin in range(Ntomo):
    cata1 = cata[cata['tom_bin_id']==(i_ZBbin+1)]
    print('+++ ZB (min, max, mean, median, Nobj, Neff)', i_ZBbin, 
          np.min(cata1['phz_mode_1']), np.max(cata1['phz_mode_1']),
          np.mean(cata1['phz_mode_1']), np.median(cata1['phz_mode_1']),
          len(cata1),
          (np.sum(cata1['shape_weight']))**2 / np.sum(cata1['shape_weight']**2)
          )
    for j_ZBbin in range(i_ZBbin, Ntomo):
        ## where to save
        outpath = os.path.join(outdir, 
                               f'XI_{shear_type}_whole_nbins_{theta_nbins}_theta_{theta_min:.1f}_{theta_max:.1f}_zbins_{i_ZBbin+1}_{j_ZBbin+1}.asc')
        ## cross correlation
        if i_ZBbin != j_ZBbin:
            cata2 = cata[cata['tom_bin_id']==(j_ZBbin+1)]
            print('+++ ZB (min, max, mean, median, Nobj, Neff)', j_ZBbin, 
                np.min(cata2['phz_mode_1']), np.max(cata2['phz_mode_1']),
                np.mean(cata2['phz_mode_1']), np.median(cata2['phz_mode_1']),
                len(cata2),
                (np.sum(cata2['shape_weight']))**2 / np.sum(cata2['shape_weight']**2)
                )
        ## auto correlation
        else:
            cata2 = None
        ## run
        CalCorrFunc(outpath,
                cata1, cata2=cata2, 
                theta_Nbin=theta_nbins, theta_bin_slop=theta_bin_slop,
                theta_min=theta_min, theta_max=theta_max, theta_unit=theta_unit, 
                theta_bin_type='Log', 
                num_threads=nthr, 
                c12_cata1=[0, 0], 
                c12_cata2=[0, 0],
                var_method='shot',
                Npatch=None,
                outpath_cov=None)
        ## extract used values
        cata_tmp = np.loadtxt(outpath)
        Nrows = len(cata_tmp[:, 0])
        cata_tmp = pd.DataFrame({'ito': (i_ZBbin * np.ones(Nrows)).astype(int),
                                'jto': (j_ZBbin * np.ones(Nrows)).astype(int),
                                'theta': cata_tmp[:, 1],
                                'xi_p': cata_tmp[:, 3],
                                'xi_m': cata_tmp[:, 4],
                                'xierr_p': cata_tmp[:, 7],
                                'xierr_m': cata_tmp[:, 8]
                                })
        cata_final.append(cata_tmp)
## combine and save
cata_final = pd.concat(cata_final)
outpath = os.path.join(outdir, f'XI_{shear_type}_whole_nbins_{theta_nbins}_theta_{theta_min:.1f}_{theta_max:.1f}_combined.csv')
cata_final.to_csv(outpath, index=False)
print('combined cata saved to', outpath)
outpath_2pt_list.append(outpath)

## subsamples
print('++++++ Measuring 2pt for the subsamples')
for i_half in range(2):
    cata_selec = cata[(cata['subsample_ID'] == i_half)].reset_index(drop=True)
    print(f"++++++ Number in subsample {i_half}", len(cata_selec))
    ## save info
    if i_half == 0:
        save_label = 'blue_lt'
    else:
        save_label = 'red_geq'
    ## loop over tomo bins
    cata_final = []
    for i_ZBbin in range(Ntomo):
        cata1 = cata_selec[cata_selec['tom_bin_id']==(i_ZBbin+1)]
        print('+++ ZB (min, max, mean, median, Nobj, Neff)', i_ZBbin, 
            np.min(cata1['phz_mode_1']), np.max(cata1['phz_mode_1']),
            np.mean(cata1['phz_mode_1']), np.median(cata1['phz_mode_1']),
            len(cata1),
            (np.sum(cata1['shape_weight']))**2 / np.sum(cata1['shape_weight']**2)
            )
        for j_ZBbin in range(i_ZBbin, Ntomo):
            ## where to save
            outpath = os.path.join(outdir,
                                f'XI_{shear_type}_U_R_{save_label}{cut_value:.1f}_nbins_{theta_nbins}_theta_{theta_min:.1f}_{theta_max:.1f}_zbins_{i_ZBbin+1}_{j_ZBbin+1}.asc')

            ## cross correlation
            if i_ZBbin != j_ZBbin:
                cata2 = cata_selec[cata_selec['tom_bin_id']==(j_ZBbin+1)]
                print('+++ ZB (min, max, mean, median, Nobj, Neff)', j_ZBbin, 
                    np.min(cata2['phz_mode_1']), np.max(cata2['phz_mode_1']),
                    np.mean(cata2['phz_mode_1']), np.median(cata2['phz_mode_1']),
                    len(cata2),
                    (np.sum(cata2['shape_weight']))**2 / np.sum(cata2['shape_weight']**2)
                    )
            ## auto correlation
            else:
                cata2 = None
            ## run
            CalCorrFunc(outpath,
                    cata1, cata2=cata2, 
                    theta_Nbin=theta_nbins, theta_bin_slop=theta_bin_slop,
                    theta_min=theta_min, theta_max=theta_max, theta_unit=theta_unit, 
                    theta_bin_type='Log', 
                    num_threads=nthr, 
                    c12_cata1=[0, 0], 
                    c12_cata2=[0, 0],
                    var_method='shot',
                    Npatch=None,
                    outpath_cov=None)
            ## extract used values
            cata_tmp = np.loadtxt(outpath)
            Nrows = len(cata_tmp[:, 0])
            cata_tmp = pd.DataFrame({'ito': (i_ZBbin * np.ones(Nrows)).astype(int),
                                    'jto': (j_ZBbin * np.ones(Nrows)).astype(int),
                                    'theta': cata_tmp[:, 1],
                                    'xi_p': cata_tmp[:, 3],
                                    'xi_m': cata_tmp[:, 4],
                                    'xierr_p': cata_tmp[:, 7],
                                    'xierr_m': cata_tmp[:, 8]
                                    })
            cata_final.append(cata_tmp)
    ## combine and save
    cata_final = pd.concat(cata_final)
    outpath = os.path.join(outdir, 
                        f'XI_{shear_type}_U_R_{save_label}{cut_value:.1f}_nbins_{theta_nbins}_theta_{theta_min:.1f}_{theta_max:.1f}_combined.csv')
    cata_final.to_csv(outpath, index=False)
    print('combined cata saved to', outpath)
    outpath_2pt_list.append(outpath)

## +++ Plot correlation function

## plot individual
CATAs = [pd.read_csv(inpath) for inpath in outpath_2pt_list]
COLORs = ['k', 'b', 'r']
LABELs = ['whole', 'U_R<1.2', 'U_R>=1.2']
yscaling = 1e4
timesX = True
ylog = False
XLABEL = r'$\theta$ [arcmin]'
XRANGE = [0.9, 410]
xlog = True
for which_half in ['xip', 'xim']:
    if which_half == 'xip':
        YLABEL = r'$\theta \times \Delta\xi_{+}$ [1e-4 arcmin]'
    else:
        YLABEL = r'$\theta \times \Delta\xi_{-}$ [1e-4 arcmin]'
    outpath = f'./plots/C_red_blue_{shear_type}_{which_half}_individual.png'
    CorrPlotFunc_half(outpath, Ntomo,
                        CATAs, which_half,
                        COLORs, 
                        yscaling = yscaling, 
                        LABELs=LABELs, 
                        LINEs=None, LINEWs=None, 
                        POINTs=None, POINTSs=None, ERRORSIZEs=None,
                        subLABEL_locX=0.1, subLABEL_locY=0.8,
                        XRANGE=XRANGE, YRANGE=None,
                        XLABEL=XLABEL, YLABEL=YLABEL, TITLE=TITLE,
                        xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                        vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                        hlines=None, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                        xlog=xlog, invertX=False, ylog=ylog, invertY=False, 
                        loc_legend='best', legend_frame=False,
                        font_size=9, usetex=False,
                        fill_between_xs_list=None, 
                        fill_between_yLows_list=None, fill_between_yHighs_list=None,
                        fill_between_COLORs_list=None, fill_between_alphas_list=None,
                        FIGSIZE=[6.4, 4.8],
                        TIGHT=False, 
                        N_rows=None, N_cols=None,
                        timesX=timesX,
                        xi_p_norm=1, xi_m_norm=1,
                        xi_p_err_col='xierr_p', xi_m_err_col='xierr_m')

## plot diff with respective to full
CATAs = [pd.read_csv(inpath) for inpath in outpath_2pt_list[1:]]
CATA_base = pd.read_csv(outpath_2pt_list[0])
COLORs = ['b', 'r']
LABELs = ['U_R<1.2', 'U_R>=1.2']
YRANGE = [-2.2, 2.2]
yscaling = 1e4
timesX = True
ylog = False
XRANGE = [0.9, 410]
XLABEL = r'$\theta$ [arcmin]'
xlog = True
hlines = [0.0]
for which_half in ['xip', 'xim']:
    if which_half == 'xip':
        YLABEL = r'$\theta \times \Delta\xi_{+}$ [1e-4 arcmin]'
    else:
        YLABEL = r'$\theta \times \Delta\xi_{-}$ [1e-4 arcmin]'
    outpath = f'./plots/C_red_blue_{shear_type}_{which_half}_diff2base.png'
    CorrPlotFunc_half_diff2base(outpath, Ntomo,
                        CATAs, CATA_base, which_half,
                        COLORs, 
                        yscaling = yscaling, 
                        LABELs=LABELs, 
                        LINEs=None, LINEWs=None, 
                        POINTs=None, POINTSs=None, ERRORSIZEs=None,
                        subLABEL_locX=0.1, subLABEL_locY=0.8,
                        XRANGE=XRANGE, YRANGE=YRANGE,
                        XLABEL=XLABEL, YLABEL=YLABEL, TITLE=TITLE,
                        xtick_min_label=True, xtick_spe=None, ytick_min_label=True, ytick_spe=None,
                        vlines=None, vline_styles=None, vline_colors=None, vline_labels=None, vline_widths=None,
                        hlines=hlines, hline_styles=None, hline_colors=None, hline_labels=None, hline_widths=None,
                        xlog=xlog, invertX=False, ylog=ylog, invertY=False, 
                        loc_legend='best', legend_frame=False,
                        font_size=9, usetex=False,
                        fill_between_xs_list=None, 
                        fill_between_yLows_list=None, fill_between_yHighs_list=None,
                        fill_between_COLORs_list=None, fill_between_alphas_list=None,
                        FIGSIZE=[6.4, 4.8],
                        TIGHT=False, 
                        N_rows=None, N_cols=None,
                        timesX=timesX,
                        xi_p_norm=1, xi_m_norm=1,
                        xi_p_err_col='xierr_p', xi_m_err_col='xierr_m')