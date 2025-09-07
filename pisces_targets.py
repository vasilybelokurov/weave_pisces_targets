#!/usr/bin/env python3
"""
WEAVE Pisces Target Selection
Converts notebook to optimized Python script for target list generation.
"""

import matplotlib as mpl
mpl.rc('font', family='serif')

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import astropy.units as u
import astropy.time
from astropy.coordinates import (SkyCoord, Distance, Galactic,
                                 EarthLocation, AltAz)
import astropy.coordinates as coord
import pandas as pd
from functools import reduce
from astropy.table import Table
import matplotlib.colors as mcolors
from matplotlib import colors
import requests
import os

# Download helper modules if not present
def download_module(url, filename):
    """Download module from GitHub if not already present."""
    if not os.path.exists(filename):
        print(f"Downloading {filename}...")
        r = requests.get(url)
        with open(filename, 'w') as f:
            f.write(r.text)
        print(f"{filename} downloaded successfully.")
    
# Download required modules
download_module('https://raw.githubusercontent.com/segasai/astrolibpy/master/astrolib/cv_coord.py', 'cv_coord.py')
download_module('https://raw.githubusercontent.com/segasai/astrolibpy/master/my_utils/sphere_rotate.py', 'sphere_rotate.py')

import cv_coord as cv_coord
import sphere_rotate as sr

# Plotting setup
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10

class PiscesTargetSelector:
    """WEAVE Pisces target selection and optimization."""
    
    def __init__(self):
        """Initialize with default survey parameters."""
        # Observing region
        self.tarx = [-45, 15]  # RA range [deg]
        self.tary = [-20, 0]   # Dec range [deg]
        
        # Target selection parameters - can be optimized
        self.rgb_params = {
            'clean': {
                'rpm_cut': 18.0,
                'par_cut': 0.2,
                'col_cut': 1.0,
                'mag_lim': 18.0,
                'dist_cut': 50,
                'priority': 10
            },
            'relaxed': {
                'rpm_cut': 18.5,
                'par_cut': 0.25,
                'col_cut': 1.0,
                'mag_lim': 18.5,
                'dist_cut': 30,
                'priority': 9
            }
        }
        
        self.bhb_params = {
            'mag_cut': 18.0,
            'priority': 8
        }
        
        self.data_loaded = False
        
    def load_catalogs(self):
        """Load all required catalogs from remote URLs."""
        print("Loading catalogs...")
        
        # RG catalog
        print("Loading RG catalog...")
        target_url = 'https://people.ast.cam.ac.uk/~vasily/data/rg_gaia2mass.fits'
        hdul = fits.open(target_url)
        self.dg = hdul[1].data
        hdul.close()
        
        # Convert RA [0,360] -> [-180,180]
        self.ra180 = self.dg['ra'].copy()
        w = np.where(self.dg['ra'] > 180)
        self.ra180[w] = self.dg['ra'][w] - 360
        
        # BHB catalog
        print("Loading BHB catalog...")
        target_url = 'https://people.ast.cam.ac.uk/~vasily/data/BHB_EDR3_2023_with_mags.fits'
        hdul = fits.open(target_url)
        self.bhb = hdul[1].data
        hdul.close()
        
        # Convert RA [0,360] -> [-180,180] for BHB
        self.bhbra180 = self.bhb['ra'].copy()
        w = np.where(self.bhb['ra'] > 180)
        self.bhbra180[w] = self.bhb['ra'][w] - 360
        
        # Gaia XP catalog (for future use)
        print("Loading Gaia XP catalog...")
        target_url = 'https://people.ast.cam.ac.uk/~vasily/data/catalogues/gaia_xp_andrae_parallax01.fits'
        hdul = fits.open(target_url)
        self.da = hdul[1].data
        hdul.close()
        
        self.data_loaded = True
        print("All catalogs loaded successfully.")
        
    def apply_global_rgb_selection(self, params_key='clean'):
        """Apply RGB star selection cuts globally (no field constraint)."""
        params = self.rgb_params[params_key]
        
        # Apply selection cuts globally
        wpar = np.where(self.dg['parallax'] < params['par_cut'])
        wcol = np.where(self.dg['col'] > params['col_cut'])
        wmag = np.where(self.dg['mg'] < params['mag_lim'])
        wdist = np.where(self.dg['dist'] > params['dist_cut'])
        wrpm = np.where(self.dg['rpm'] < params['rpm_cut'])
        
        # Combine all cuts
        wtar_global = reduce(np.intersect1d, (wpar, wcol, wmag, wdist, wrpm))
        
        return wtar_global
        
    def apply_global_bhb_selection(self):
        """Apply BHB star selection cuts globally (no field constraint)."""
        # Distance cut based on magnitude
        wbhb_dist = np.where(self.bhb['ps1g'] > self.bhb_params['mag_cut'])
        
        return wbhb_dist[0]

    def apply_rgb_selection(self, params_key='clean'):
        """Apply RGB star selection cuts."""
        params = self.rgb_params[params_key]
        
        # Define observing region
        wxy = np.where((self.ra180 > self.tarx[0]) & 
                      (self.ra180 < self.tarx[1]) & 
                      (self.dg['dec'] > self.tary[0]) & 
                      (self.dg['dec'] < self.tary[1]))
        
        # Apply selection cuts
        wpar = np.where(self.dg['parallax'] < params['par_cut'])
        wcol = np.where(self.dg['col'] > params['col_cut'])
        wmag = np.where(self.dg['mg'] < params['mag_lim'])
        wdist = np.where(self.dg['dist'] > params['dist_cut'])
        wrpm = np.where(self.dg['rpm'] < params['rpm_cut'])
        
        # Combine all cuts
        wtar = reduce(np.intersect1d, (wpar, wcol, wmag, wrpm, wdist))
        wtar_xy = np.intersect1d(wtar, wxy[0])
        
        return wtar_xy
        
    def apply_bhb_selection(self):
        """Apply BHB star selection cuts."""
        # Define observing region
        wbhb_xy = np.where((self.bhbra180 > self.tarx[0]) & 
                          (self.bhbra180 < self.tarx[1]) & 
                          (self.bhb['dec'] > self.tary[0]) & 
                          (self.bhb['dec'] < self.tary[1]))
        
        # Distance cut based on magnitude
        wbhb_dist = np.where(self.bhb['ps1g'] > self.bhb_params['mag_cut'])
        
        wtar2_xy = np.intersect1d(wbhb_dist[0], wbhb_xy[0])
        
        return wtar2_xy
        
    def generate_target_list(self, date_suffix=None):
        """Generate complete target list with priorities."""
        if not self.data_loaded:
            self.load_catalogs()
            
        print("\nGenerating target list...")
        
        # Get clean RGB targets (highest priority)
        wtar1_xy = self.apply_rgb_selection('clean')
        print(f'Clean RGB targets: {len(wtar1_xy)}')
        
        # Get relaxed RGB targets
        wtar12_xy = self.apply_rgb_selection('relaxed')
        print(f'Total relaxed RGB targets: {len(wtar12_xy)}')
        
        # Remove overlap with clean targets
        keep_wtar12_xy = list(set(wtar12_xy) - set(wtar1_xy))
        keep_wtar12_xy = np.array(keep_wtar12_xy)
        print(f'Additional relaxed RGB targets: {len(keep_wtar12_xy)}')
        
        # Get BHB targets
        wtar2_xy = self.apply_bhb_selection()
        print(f'BHB targets: {len(wtar2_xy)}')
        
        # Create target arrays
        ntar1 = len(wtar1_xy)
        ntar12 = len(keep_wtar12_xy)
        ntar2 = len(wtar2_xy)
        
        # Pack Target 1 (Clean RGB)
        rev_cur1 = np.full(ntar1, 3, dtype=np.int64)
        ps1_cur1 = np.full(ntar1, 0, dtype=np.int64)
        pri_cur1 = np.full(ntar1, self.rgb_params['clean']['priority'], dtype=np.int64)
        dtar1 = np.rec.fromarrays([
            self.dg['source_id'][wtar1_xy],
            self.dg['ra'][wtar1_xy], 
            self.dg['dec'][wtar1_xy], 
            rev_cur1, ps1_cur1, pri_cur1
        ], names=['SOURCE_ID', 'RA', 'DEC', 'GAIA_REV_ID', 'PS1_ID', 'PRIORITY'])
        
        # Pack Target 12 (Relaxed RGB)
        rev_cur12 = np.full(ntar12, 3, dtype=np.int64)
        ps1_cur12 = np.full(ntar12, 0, dtype=np.int64)
        pri_cur12 = np.full(ntar12, self.rgb_params['relaxed']['priority'], dtype=np.int64)
        dtar12 = np.rec.fromarrays([
            self.dg['source_id'][keep_wtar12_xy],
            self.dg['ra'][keep_wtar12_xy], 
            self.dg['dec'][keep_wtar12_xy], 
            rev_cur12, ps1_cur12, pri_cur12
        ], names=['SOURCE_ID', 'RA', 'DEC', 'GAIA_REV_ID', 'PS1_ID', 'PRIORITY'])
        
        # Pack Target 2 (BHB)
        rev_cur2 = np.full(ntar2, 3, dtype=np.int64)
        ps1_cur2 = np.full(ntar2, 0, dtype=np.int64)
        pri_cur2 = np.full(ntar2, self.bhb_params['priority'], dtype=np.int64)
        dtar2 = np.rec.fromarrays([
            self.bhb['source_id'][wtar2_xy],
            self.bhb['ra'][wtar2_xy], 
            self.bhb['dec'][wtar2_xy], 
            rev_cur2, ps1_cur2, pri_cur2
        ], names=['SOURCE_ID', 'RA', 'DEC', 'GAIA_REV_ID', 'PS1_ID', 'PRIORITY'])
        
        # Concatenate all targets
        data = np.concatenate((dtar1, dtar12, dtar2))
        
        print(f'\nTotal targets: {len(data)}')
        print(f'  - Clean RGB (priority {self.rgb_params["clean"]["priority"]}): {ntar1}')
        print(f'  - Relaxed RGB (priority {self.rgb_params["relaxed"]["priority"]}): {ntar12}')
        print(f'  - BHB (priority {self.bhb_params["priority"]}): {ntar2}')
        
        return data, (wtar1_xy, keep_wtar12_xy, wtar2_xy)
        
    def save_target_list(self, data, date_suffix=None):
        """Save target list to FITS file."""
        if date_suffix is None:
            date_suffix = '021023'  # Default from notebook
            
        filename = f'Pisces_{date_suffix}.fits'
        t = Table(data, copy=False)
        t.write(filename, overwrite=True)
        
        # Add metadata
        with fits.open(filename, mode='update') as hdul:
            hdul[0].header['VERSION'] = date_suffix
            hdul[0].header['SURVEY'] = 'WEAVE Pisces'
            hdul[0].header['TARGET_N'] = len(data)
            
        print(f'Target list saved to: {filename}')
        return filename
        
    def plot_all_diagnostics(self, indices_list, show_plots=True):
        """Create all diagnostic plots from the original notebook."""
        wtar1_xy, keep_wtar12_xy, wtar2_xy = indices_list
        
        # Get global selections (before field constraint) - this was my mistake!
        wtar1_global = self.apply_global_rgb_selection('clean')
        wtar12_global = self.apply_global_rgb_selection('relaxed') 
        wtar2_global = self.apply_global_bhb_selection()
        
        # 1. Full sky distribution of clean RGB targets (ALL targets globally)
        plt.figure(figsize=[10,5])
        plt.plot(self.ra180[wtar1_global], self.dg['dec'][wtar1_global], 'o', 
                color='black', markersize=4, alpha=0.5, label=f'Clean RGB ({len(wtar1_global)} global)')
        plt.plot([self.tarx[0], self.tarx[1], self.tarx[1], self.tarx[0], self.tarx[0]],
                [self.tary[0], self.tary[0], self.tary[1], self.tary[1], self.tary[0]],
                color='red', linestyle='solid', linewidth=2, label='Survey region')
        plt.xlim(180, -180)
        plt.ylim(-90, 90)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'Clean RGB Targets - Full Sky (Global: {len(wtar1_global)}, In field: {len(wtar1_xy)})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('plots/pisces_clean_rgb_fullsky.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 2. Full sky distribution of relaxed RGB targets (ALL targets globally)
        plt.figure(figsize=[10,5])
        plt.plot(self.ra180[wtar12_global], self.dg['dec'][wtar12_global], 'o', 
                color='black', markersize=3, alpha=0.1, label=f'Relaxed RGB ({len(wtar12_global)} global)')
        plt.plot([self.tarx[0], self.tarx[1], self.tarx[1], self.tarx[0], self.tarx[0]],
                [self.tary[0], self.tary[0], self.tary[1], self.tary[1], self.tary[0]],
                color='red', linestyle='solid', linewidth=2, label='Survey region')
        plt.xlim(180, -180)
        plt.ylim(-90, 90)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'Relaxed RGB Targets - Full Sky (Global: {len(wtar12_global)}, In field: {len(wtar1_xy)+len(keep_wtar12_xy)})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('plots/pisces_relaxed_rgb_fullsky.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 3. Sagittarius coordinate system plot (using global clean RGB)
        phi1, phi2 = sr.sphere_rotate(self.dg['ra'], self.dg['dec'], 72., -14., 191.10487)
        plt.figure(figsize=[10,5])
        plt.plot(phi1[wtar1_global], phi2[wtar1_global], 'o', 
                color='black', markersize=4, alpha=0.3, label=f'Clean RGB ({len(wtar1_global)})')
        plt.xlim(180, -180)
        plt.ylim(-40, 40)
        plt.xlabel('Phi1 [deg] (Sgr coords)')
        plt.ylabel('Phi2 [deg] (Sgr coords)')
        plt.title('Clean RGB in Sagittarius Coordinates')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('pisces_clean_rgb_sgr_coords.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 4. 2D histogram in Sagittarius coordinates (using global clean RGB)
        plt.figure(figsize=[10,5])
        plt.hist2d(phi1[wtar1_global], phi2[wtar1_global], bins=[50,100], vmax=30)
        plt.colorbar(label='Target count')
        plt.ylim(-20, 20)
        plt.xlabel('Phi1 [deg] (Sgr coords)')
        plt.ylabel('Phi2 [deg] (Sgr coords)')
        plt.title('Clean RGB Density in Sagittarius Coordinates')
        if show_plots:
            plt.savefig('plots/pisces_clean_rgb_sgr_density.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 5. BHB targets full sky (ALL targets globally)
        plt.figure(figsize=[10,5])
        plt.plot(self.bhbra180[wtar2_global], self.bhb['dec'][wtar2_global], 'o', 
                color='black', markersize=2, alpha=0.25, label=f'BHB ({len(wtar2_global)} global)')
        plt.plot([self.tarx[0], self.tarx[1], self.tarx[1], self.tarx[0], self.tarx[0]],
                [self.tary[0], self.tary[0], self.tary[1], self.tary[1], self.tary[0]],
                color='red', linestyle='solid', linewidth=2, label='Survey region')
        plt.xlim(180, -180)
        plt.ylim(-90, 90)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'BHB Targets - Full Sky (Global: {len(wtar2_global)}, In field: {len(wtar2_xy)})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('pisces_bhb_fullsky.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 6. Survey region zoom - Clean RGB
        plt.figure(figsize=[10,5])
        plt.plot(self.ra180[wtar1_xy], self.dg['dec'][wtar1_xy], 'o', 
                color='black', markersize=4, alpha=0.5, label='Clean RGB')
        plt.xlim(self.tarx)
        plt.ylim(self.tary)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'Clean RGB - Survey Region Detail ({len(wtar1_xy)} targets)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('plots/pisces_clean_rgb_zoom.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 7. Survey region zoom - Relaxed RGB (additional only)
        plt.figure(figsize=[10,5])
        plt.plot(self.ra180[keep_wtar12_xy], self.dg['dec'][keep_wtar12_xy], 'o', 
                color='blue', markersize=3, alpha=0.5, label='Additional relaxed RGB')
        plt.xlim(self.tarx)
        plt.ylim(self.tary)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'Additional Relaxed RGB - Survey Region ({len(keep_wtar12_xy)} targets)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('pisces_relaxed_rgb_zoom.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 8. Survey region zoom - BHB
        plt.figure(figsize=[10,5])
        plt.plot(self.bhbra180[wtar2_xy], self.bhb['dec'][wtar2_xy], 'o', 
                color='green', markersize=3, alpha=0.5, label='BHB targets')
        plt.xlim(self.tarx)
        plt.ylim(self.tary)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'BHB Targets - Survey Region ({len(wtar2_xy)} targets)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('plots/pisces_bhb_zoom.png', dpi=150, bbox_inches='tight')
            plt.show()
        
        # 9. Combined survey region plot
        plt.figure(figsize=[12,8])
        plt.plot(self.ra180[wtar1_xy], self.dg['dec'][wtar1_xy], 'o', 
                color='red', markersize=4, alpha=0.7, label=f'Clean RGB ({len(wtar1_xy)})')
        plt.plot(self.ra180[keep_wtar12_xy], self.dg['dec'][keep_wtar12_xy], 'o', 
                color='blue', markersize=3, alpha=0.5, label=f'Additional relaxed RGB ({len(keep_wtar12_xy)})')
        plt.plot(self.bhbra180[wtar2_xy], self.bhb['dec'][wtar2_xy], 'o', 
                color='green', markersize=3, alpha=0.5, label=f'BHB ({len(wtar2_xy)})')
        plt.xlim(self.tarx)
        plt.ylim(self.tary)
        plt.xlabel('RA [deg]')
        plt.ylabel('Dec [deg]')
        plt.title(f'All Targets - Survey Region (Total: {len(wtar1_xy) + len(keep_wtar12_xy) + len(wtar2_xy)})')
        plt.legend()
        plt.grid(True, alpha=0.3)
        if show_plots:
            plt.savefig('plots/pisces_all_targets_combined.png', dpi=150, bbox_inches='tight')
            plt.show()
            
        print(f"\nDiagnostic plots saved:")
        print("- pisces_clean_rgb_fullsky.png")
        print("- pisces_relaxed_rgb_fullsky.png") 
        print("- pisces_clean_rgb_sgr_coords.png")
        print("- pisces_clean_rgb_sgr_density.png")
        print("- pisces_bhb_fullsky.png")
        print("- pisces_clean_rgb_zoom.png")
        print("- pisces_relaxed_rgb_zoom.png")
        print("- pisces_bhb_zoom.png")
        print("- pisces_all_targets_combined.png")
        
    def plot_targets(self, indices_list, labels=None, colors=None, save_fig=None):
        """Simple target distribution plot - kept for compatibility."""
        return self.plot_all_diagnostics(indices_list, show_plots=True)
        
    def generate_target_summary_table(self, n_pointings=35):
        """Generate comprehensive target summary table."""
        print("\n" + "="*80)
        print("WEAVE PISCES TARGET SELECTION SUMMARY")
        print("="*80)
        
        # Get global and field-constrained selections
        wtar1_global = self.apply_global_rgb_selection('clean')
        wtar12_global = self.apply_global_rgb_selection('relaxed')
        wtar2_global = self.apply_global_bhb_selection()
        
        wtar1_xy = self.apply_rgb_selection('clean')
        wtar12_xy = self.apply_rgb_selection('relaxed')
        wtar2_xy = self.apply_bhb_selection()
        
        # Remove overlap between clean and relaxed RGB
        keep_wtar12_xy = list(set(wtar12_xy) - set(wtar1_xy))
        keep_wtar12_xy = np.array(keep_wtar12_xy)
        
        # Survey area calculations
        survey_area_deg2 = (self.tarx[1] - self.tarx[0]) * (self.tary[1] - self.tary[0])
        pointing_area_deg2 = np.pi * (2.0/2)**2  # WEAVE 2° diameter field
        
        print(f"Survey Region: RA=[{self.tarx[0]}, {self.tarx[1]}]°, Dec=[{self.tary[0]}, {self.tary[1]}]°")
        print(f"Survey Area: {survey_area_deg2:.1f} deg²")
        print(f"WEAVE Field of View: {pointing_area_deg2:.2f} deg² (2° diameter)")
        print(f"Requested Pointings: {n_pointings}")
        print(f"Total Coverage: {n_pointings * pointing_area_deg2:.1f} deg²")
        
        # Create summary table
        print("\n" + "-"*120)
        print(f"{'TARGET CATEGORY':<25} {'PRIORITY':<8} {'GLOBAL':<8} {'IN FIELD':<9} {'PER FIELD':<10} {'SELECTION CRITERIA':<50}")
        print("-"*120)
        
        # Clean RGB
        criteria_clean = (f"parallax<{self.rgb_params['clean']['par_cut']}, "
                         f"RPM<{self.rgb_params['clean']['rpm_cut']}, "
                         f"(BP-RP)₀>{self.rgb_params['clean']['col_cut']}, "
                         f"G₀<{self.rgb_params['clean']['mag_lim']}, "
                         f"dist>{self.rgb_params['clean']['dist_cut']}kpc")
        
        print(f"{'Clean RGB':<25} {self.rgb_params['clean']['priority']:<8} {len(wtar1_global):<8} "
              f"{len(wtar1_xy):<9} {len(wtar1_xy)/n_pointings:<10.1f} {criteria_clean:<50}")
        
        # Relaxed RGB (additional only)
        criteria_relaxed = (f"parallax<{self.rgb_params['relaxed']['par_cut']}, "
                           f"RPM<{self.rgb_params['relaxed']['rpm_cut']}, "
                           f"(BP-RP)₀>{self.rgb_params['relaxed']['col_cut']}, "
                           f"G₀<{self.rgb_params['relaxed']['mag_lim']}, "
                           f"dist>{self.rgb_params['relaxed']['dist_cut']}kpc")
        
        print(f"{'Additional Relaxed RGB':<25} {self.rgb_params['relaxed']['priority']:<8} "
              f"{len(wtar12_global)-len(wtar1_global):<8} {len(keep_wtar12_xy):<9} "
              f"{len(keep_wtar12_xy)/n_pointings:<10.1f} {criteria_relaxed:<50}")
        
        # BHB
        criteria_bhb = f"PS1_g>{self.bhb_params['mag_cut']} (distance tracer)"
        
        print(f"{'Blue Horizontal Branch':<25} {self.bhb_params['priority']:<8} {len(wtar2_global):<8} "
              f"{len(wtar2_xy):<9} {len(wtar2_xy)/n_pointings:<10.1f} {criteria_bhb:<50}")
        
        print("-"*120)
        
        # Totals
        total_global = len(wtar1_global) + (len(wtar12_global)-len(wtar1_global)) + len(wtar2_global)
        total_field = len(wtar1_xy) + len(keep_wtar12_xy) + len(wtar2_xy)
        
        print(f"{'TOTAL':<25} {'-':<8} {total_global:<8} {total_field:<9} "
              f"{total_field/n_pointings:<10.1f} {'Combined selection':<50}")
        
        print("-"*120)
        
        # Target density analysis
        print(f"\nTARGET DENSITY ANALYSIS:")
        print(f"Current density: {total_field/n_pointings:.1f} targets per pointing")
        print(f"Target density:  40.0 targets per pointing")
        print(f"Density ratio:   {(total_field/n_pointings)/40.0:.2f}× higher than target")
        
        # Efficiency metrics
        print(f"\nSELECTION EFFICIENCY:")
        print(f"Clean RGB:       {len(wtar1_xy)/len(wtar1_global)*100:.1f}% of global sample in field")
        print(f"Relaxed RGB:     {len(wtar12_xy)/len(wtar12_global)*100:.1f}% of global sample in field")
        print(f"BHB:             {len(wtar2_xy)/len(wtar2_global)*100:.1f}% of global sample in field")
        print(f"Overall:         {total_field/total_global*100:.1f}% of global targets in field")
        
        # Create detailed breakdown by magnitude and distance
        print(f"\nTARGET PROPERTIES (in survey field only):")
        print(f"Clean RGB magnitude range:     {self.dg['mg'][wtar1_xy].min():.1f} < G₀ < {self.dg['mg'][wtar1_xy].max():.1f}")
        print(f"Clean RGB distance range:      {self.dg['dist'][wtar1_xy].min():.0f} < d < {self.dg['dist'][wtar1_xy].max():.0f} kpc")
        print(f"Relaxed RGB magnitude range:   {self.dg['mg'][keep_wtar12_xy].min():.1f} < G₀ < {self.dg['mg'][keep_wtar12_xy].max():.1f}")
        print(f"Relaxed RGB distance range:    {self.dg['dist'][keep_wtar12_xy].min():.0f} < d < {self.dg['dist'][keep_wtar12_xy].max():.0f} kpc")
        print(f"BHB PS1 g range:               {self.bhb['ps1g'][wtar2_xy].min():.1f} < g < {self.bhb['ps1g'][wtar2_xy].max():.1f}")
        
        return {
            'global_counts': {'clean_rgb': len(wtar1_global), 
                            'relaxed_rgb': len(wtar12_global), 
                            'bhb': len(wtar2_global)},
            'field_counts': {'clean_rgb': len(wtar1_xy),
                           'additional_relaxed_rgb': len(keep_wtar12_xy),
                           'bhb': len(wtar2_xy)},
            'total_field': total_field,
            'density_per_pointing': total_field/n_pointings,
            'survey_area': survey_area_deg2,
            'efficiency': total_field/total_global*100
        }

    def optimize_selection(self, target_per_pointing=40, n_pointings=35):
        """Optimize target selection parameters to achieve desired target density."""
        total_target = target_per_pointing * n_pointings
        print(f"Target: {total_target} targets ({target_per_pointing} per pointing, {n_pointings} pointings)")
        
        # Generate current target list
        data, indices_list = self.generate_target_list()
        current_total = len(data)
        
        # Generate summary table
        summary_stats = self.generate_target_summary_table(n_pointings)
        
        return data, indices_list, summary_stats

def main():
    """Main execution function."""
    print("WEAVE Pisces Target Selection")
    print("=" * 40)
    
    # Initialize selector
    selector = PiscesTargetSelector()
    
    # Generate targets with summary
    data, indices_list, summary_stats = selector.optimize_selection()
    
    # Create plots
    selector.plot_targets(indices_list, save_fig='plots/pisces_targets_distribution.png')
    
    # Save target list
    filename = selector.save_target_list(data)
    
    print(f"\nTarget selection complete!")
    print(f"Output file: {filename}")

if __name__ == "__main__":
    main()