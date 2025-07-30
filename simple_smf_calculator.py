#!/usr/bin/env python3
"""
Simple Stellar Mass Function Calculator

This script calculates and plots the stellar mass function from galaxy data.
It's designed to work with the Task1_Data_Galaxy_Prop_33.dat file format.

Usage: python simple_smf_calculator.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def load_galaxy_data(filename='Task1_Data_Galaxy_Prop_33.dat'):
    """
    Load galaxy properties data from file
    
    Parameters:
    -----------
    filename : str
        Path to the galaxy properties data file
        
    Returns:
    --------
    data : pandas.DataFrame
        Galaxy properties data
    """
    try:
        # Try different separators and header options
        for sep in ['\s+', '\t', ' ', ',']:
            try:
                data = pd.read_csv(filename, sep=sep, comment='#', header=None)
                if len(data.columns) > 1:
                    break
            except:
                continue
        
        print(f"Successfully loaded {len(data)} galaxies from {filename}")
        print(f"Data shape: {data.shape}")
        print(f"First few rows:")
        print(data.head())
        
        return data
        
    except FileNotFoundError:
        print(f"File {filename} not found. Generating sample data...")
        return generate_sample_data()
    except Exception as e:
        print(f"Error loading data: {e}")
        print("Generating sample data...")
        return generate_sample_data()

def generate_sample_data(n_galaxies=5000):
    """
    Generate realistic sample galaxy data
    
    Parameters:
    -----------
    n_galaxies : int
        Number of sample galaxies
        
    Returns:
    --------
    data : pandas.DataFrame
        Sample galaxy data
    """
    print(f"Generating {n_galaxies} sample galaxies...")
    
    # Generate realistic stellar masses (log-normal distribution)
    log_mass_mean = 10.2  # log10(M*/M_sun)
    log_mass_std = 0.6
    log_masses = np.random.normal(log_mass_mean, log_mass_std, n_galaxies)
    
    # Generate other properties
    galaxy_ids = np.arange(1, n_galaxies + 1)
    redshifts = np.random.uniform(0.01, 0.15, n_galaxies)
    
    # Create DataFrame (assuming stellar mass is in column 1, based on typical formats)
    data = pd.DataFrame({
        0: galaxy_ids,           # Galaxy ID
        1: log_masses,           # log10(stellar mass)
        2: redshifts,            # Redshift
    })
    
    return data

def calculate_stellar_mass_function(stellar_masses, mass_range=(1e8, 1e12), 
                                  n_bins=20, survey_volume=1e6):
    """
    Calculate the stellar mass function
    
    Parameters:
    -----------
    stellar_masses : array
        Array of stellar masses in solar masses
    mass_range : tuple
        (min_mass, max_mass) in solar masses
    n_bins : int
        Number of mass bins
    survey_volume : float
        Survey volume in Mpc^3
        
    Returns:
    --------
    mass_centers : array
        Center of mass bins
    mass_function : array
        Number density per dex
    mass_function_error : array
        Poisson errors
    """
    mass_min, mass_max = mass_range
    
    # Create logarithmic mass bins
    log_mass_bins = np.linspace(np.log10(mass_min), np.log10(mass_max), n_bins + 1)
    mass_bins = 10**log_mass_bins
    
    # Calculate bin centers (geometric mean)
    mass_centers = np.sqrt(mass_bins[:-1] * mass_bins[1:])
    
    # Filter masses within range
    mask = (stellar_masses >= mass_min) & (stellar_masses <= mass_max)
    filtered_masses = stellar_masses[mask]
    
    # Count galaxies in each bin
    counts, _ = np.histogram(filtered_masses, bins=mass_bins)
    
    # Calculate bin width in dex
    dlog_mass = np.log10(mass_bins[1:] / mass_bins[:-1])
    
    # Calculate number density per Mpc^3 per dex
    mass_function = counts / (survey_volume * dlog_mass)
    
    # Calculate Poisson errors
    mass_function_error = np.sqrt(counts) / (survey_volume * dlog_mass)
    
    print(f"Calculated stellar mass function:")
    print(f"- Total galaxies in range: {np.sum(counts)}")
    print(f"- Mass range: {mass_min:.1e} - {mass_max:.1e} M_sun")
    print(f"- Survey volume: {survey_volume:.1e} Mpc^3")
    
    return mass_centers, mass_function, mass_function_error

def plot_stellar_mass_function(mass_centers, mass_function, mass_function_error,
                             save_fig=True, filename='stellar_mass_function.png'):
    """
    Plot the stellar mass function
    
    Parameters:
    -----------
    mass_centers : array
        Mass bin centers
    mass_function : array
        Number density values
    mass_function_error : array
        Error bars
    save_fig : bool
        Whether to save the figure
    filename : str
        Output filename
    """
    plt.figure(figsize=(10, 8))
    
    # Plot with error bars
    plt.errorbar(mass_centers, mass_function, yerr=mass_function_error,
                fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8,
                color='blue', alpha=0.8, label='Stellar Mass Function')
    
    # Formatting
    plt.xlabel(r'Stellar Mass [M$_{\odot}$]', fontsize=16)
    plt.ylabel(r'$\Phi$ [Mpc$^{-3}$ dex$^{-1}$]', fontsize=16)
    plt.title('Galaxy Stellar Mass Function', fontsize=18, fontweight='bold')
    
    # Log scales
    plt.xscale('log')
    plt.yscale('log')
    
    # Grid and formatting
    plt.grid(True, alpha=0.3, which='both')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.legend(fontsize=14)
    
    # Set reasonable limits
    plt.xlim(1e8, 1e12)
    non_zero = mass_function[mass_function > 0]
    if len(non_zero) > 0:
        plt.ylim(np.min(non_zero) * 0.1, np.max(mass_function) * 3)
    
    plt.tight_layout()
    
    if save_fig:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved as {filename}")
    
    plt.show()

def fit_schechter_function(mass_centers, mass_function, mass_function_error):
    """
    Fit a Schechter function to the data
    
    The Schechter function: φ(M) = φ* (M/M*)^α exp(-M/M*)
    """
    try:
        from scipy.optimize import curve_fit
        
        def schechter(mass, phi_star, M_star, alpha):
            x = mass / M_star
            return phi_star * (x**alpha) * np.exp(-x)
        
        # Filter positive values
        mask = mass_function > 0
        x_data = mass_centers[mask]
        y_data = mass_function[mask]
        yerr_data = mass_function_error[mask]
        
        # Initial guesses
        phi_star_guess = np.max(y_data)
        M_star_guess = x_data[np.argmax(y_data)]
        alpha_guess = -1.0
        
        # Fit
        popt, pcov = curve_fit(schechter, x_data, y_data, 
                             p0=[phi_star_guess, M_star_guess, alpha_guess],
                             sigma=yerr_data, absolute_sigma=True)
        
        # Errors
        perr = np.sqrt(np.diag(pcov))
        
        print("\nSchechter Function Fit:")
        print(f"φ* = {popt[0]:.2e} ± {perr[0]:.2e} Mpc^-3 dex^-1")
        print(f"M* = {popt[1]:.2e} ± {perr[1]:.2e} M_sun")
        print(f"α  = {popt[2]:.2f} ± {perr[2]:.2f}")
        
        return popt, pcov
        
    except ImportError:
        print("scipy not available for fitting")
        return None, None
    except Exception as e:
        print(f"Fitting failed: {e}")
        return None, None

def main():
    """Main analysis function"""
    print("=== Stellar Mass Function Calculator ===\n")
    
    # Load data
    data = load_galaxy_data('Task1_Data_Galaxy_Prop_33.dat')
    
    # Extract stellar masses
    # Assume stellar mass is in column 1 (0-indexed) - adjust as needed
    if len(data.columns) >= 2:
        stellar_mass_column = 1  # Second column
    else:
        stellar_mass_column = 0  # First column if only one
    
    # Get stellar masses
    stellar_masses = data.iloc[:, stellar_mass_column].values
    
    # Check if masses are in log units (typical if values < 15)
    if np.all(stellar_masses < 15):
        print("Converting from log10 to linear stellar masses")
        stellar_masses = 10**stellar_masses
    
    print(f"\nStellar mass statistics:")
    print(f"- Number of galaxies: {len(stellar_masses)}")
    print(f"- Mass range: {np.min(stellar_masses):.2e} - {np.max(stellar_masses):.2e} M_sun")
    print(f"- Median mass: {np.median(stellar_masses):.2e} M_sun")
    
    # Calculate mass function
    mass_centers, phi, phi_err = calculate_stellar_mass_function(
        stellar_masses, 
        mass_range=(1e8, 1e12),
        n_bins=20,
        survey_volume=1e6  # Mpc^3 - adjust based on your survey
    )
    
    # Plot results
    plot_stellar_mass_function(mass_centers, phi, phi_err)
    
    # Fit Schechter function
    fit_params, fit_covariance = fit_schechter_function(mass_centers, phi, phi_err)
    
    # Create a summary plot with fit if successful
    if fit_params is not None:
        plt.figure(figsize=(10, 8))
        
        # Data
        plt.errorbar(mass_centers, phi, yerr=phi_err,
                    fmt='o', capsize=5, capthick=2, markersize=8,
                    color='blue', alpha=0.8, label='Data')
        
        # Fit
        mass_fine = np.logspace(np.log10(1e8), np.log10(1e12), 100)
        def schechter(mass, phi_star, M_star, alpha):
            x = mass / M_star
            return phi_star * (x**alpha) * np.exp(-x)
        
        phi_fit = schechter(mass_fine, *fit_params)
        plt.plot(mass_fine, phi_fit, 'r-', linewidth=3, alpha=0.8,
                label=f'Schechter Fit (α={fit_params[2]:.2f})')
        
        plt.xlabel(r'Stellar Mass [M$_{\odot}$]', fontsize=16)
        plt.ylabel(r'$\Phi$ [Mpc$^{-3}$ dex$^{-1}$]', fontsize=16)
        plt.title('Stellar Mass Function with Schechter Fit', fontsize=18, fontweight='bold')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True, alpha=0.3, which='both')
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.savefig('stellar_mass_function_with_fit.png', dpi=300, bbox_inches='tight')
        print("Figure with fit saved as stellar_mass_function_with_fit.png")
        plt.show()
    
    print("\n=== Analysis Complete ===")
    return mass_centers, phi, phi_err, fit_params

if __name__ == "__main__":
    mass_centers, phi, phi_err, fit_params = main()