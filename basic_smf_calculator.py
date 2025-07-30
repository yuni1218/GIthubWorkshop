#!/usr/bin/env python3
"""
Basic Stellar Mass Function Calculator

This script calculates and plots the stellar mass function using only standard Python libraries.
It's designed to work without external dependencies for testing purposes.

Usage: python3 basic_smf_calculator.py
"""

import math
import random

def load_galaxy_data(filename='Task1_Data_Galaxy_Prop_33.dat'):
    """
    Load galaxy properties data from file
    
    Returns:
    --------
    data : list of lists
        Galaxy properties data
    """
    try:
        with open(filename, 'r') as f:
            data = []
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Split by whitespace and convert to float
                    row = [float(x) for x in line.split()]
                    data.append(row)
        
        print(f"Successfully loaded {len(data)} galaxies from {filename}")
        if data:
            print(f"Number of columns: {len(data[0])}")
            print(f"First few rows:")
            for i, row in enumerate(data[:5]):
                print(f"  {i}: {row}")
        
        return data
        
    except FileNotFoundError:
        print(f"File {filename} not found. Generating sample data...")
        return generate_sample_data()
    except Exception as e:
        print(f"Error loading data: {e}")
        print("Generating sample data...")
        return generate_sample_data()

def generate_sample_data(n_galaxies=1000):
    """
    Generate realistic sample galaxy data
    
    Parameters:
    -----------
    n_galaxies : int
        Number of sample galaxies
        
    Returns:
    --------
    data : list of lists
        Sample galaxy data [id, log_mass, redshift]
    """
    print(f"Generating {n_galaxies} sample galaxies...")
    
    data = []
    for i in range(n_galaxies):
        # Generate realistic stellar masses (log-normal distribution approximation)
        log_mass = random.gauss(10.2, 0.6)  # mean=10.2, std=0.6
        redshift = random.uniform(0.01, 0.15)
        
        row = [i+1, log_mass, redshift]
        data.append(row)
    
    print("Sample data generated successfully")
    return data

def calculate_stellar_mass_function(stellar_masses, mass_min=1e8, mass_max=1e12, 
                                  n_bins=20, survey_volume=1e6):
    """
    Calculate the stellar mass function
    
    Parameters:
    -----------
    stellar_masses : list
        List of stellar masses in solar masses
    mass_min : float
        Minimum stellar mass [M_sun]
    mass_max : float
        Maximum stellar mass [M_sun]
    n_bins : int
        Number of mass bins
    survey_volume : float
        Survey volume in Mpc^3
        
    Returns:
    --------
    tuple : (mass_centers, mass_function, mass_function_error)
    """
    
    # Create logarithmic mass bins
    log_mass_min = math.log10(mass_min)
    log_mass_max = math.log10(mass_max)
    
    log_mass_bins = []
    for i in range(n_bins + 1):
        log_mass = log_mass_min + (log_mass_max - log_mass_min) * i / n_bins
        log_mass_bins.append(log_mass)
    
    mass_bins = [10**log_mass for log_mass in log_mass_bins]
    
    # Calculate bin centers (geometric mean)
    mass_centers = []
    for i in range(len(mass_bins) - 1):
        center = math.sqrt(mass_bins[i] * mass_bins[i+1])
        mass_centers.append(center)
    
    # Filter masses within range
    filtered_masses = [m for m in stellar_masses if mass_min <= m <= mass_max]
    
    # Count galaxies in each bin
    counts = [0] * n_bins
    for mass in filtered_masses:
        for i in range(len(mass_bins) - 1):
            if mass_bins[i] <= mass < mass_bins[i+1]:
                counts[i] += 1
                break
    
    # Calculate bin width in dex
    dlog_mass = []
    for i in range(len(mass_bins) - 1):
        dlog = math.log10(mass_bins[i+1] / mass_bins[i])
        dlog_mass.append(dlog)
    
    # Calculate number density per Mpc^3 per dex
    mass_function = []
    mass_function_error = []
    
    for i in range(n_bins):
        phi = counts[i] / (survey_volume * dlog_mass[i])
        phi_err = math.sqrt(counts[i]) / (survey_volume * dlog_mass[i])
        
        mass_function.append(phi)
        mass_function_error.append(phi_err)
    
    print(f"Calculated stellar mass function:")
    print(f"- Total galaxies in range: {sum(counts)}")
    print(f"- Mass range: {mass_min:.1e} - {mass_max:.1e} M_sun")
    print(f"- Survey volume: {survey_volume:.1e} Mpc^3")
    
    return mass_centers, mass_function, mass_function_error

def print_results(mass_centers, mass_function, mass_function_error):
    """
    Print the results in a formatted table
    """
    print("\n" + "="*80)
    print("STELLAR MASS FUNCTION RESULTS")
    print("="*80)
    print(f"{'Mass [M_sun]':>15} {'Phi [Mpc^-3 dex^-1]':>20} {'Error':>15}")
    print("-" * 80)
    
    for i in range(len(mass_centers)):
        mass = mass_centers[i]
        phi = mass_function[i]
        err = mass_function_error[i]
        
        if phi > 0:
            print(f"{mass:>15.2e} {phi:>20.2e} {err:>15.2e}")
    
    print("-" * 80)

def create_simple_plot_data(mass_centers, mass_function):
    """
    Create simple ASCII plot data
    """
    print("\n" + "="*60)
    print("STELLAR MASS FUNCTION (ASCII PLOT)")
    print("="*60)
    
    # Filter out zero values for plotting
    plot_data = [(m, p) for m, p in zip(mass_centers, mass_function) if p > 0]
    
    if not plot_data:
        print("No data to plot")
        return
    
    # Normalize for ASCII plotting
    max_phi = max(p for m, p in plot_data)
    
    print(f"{'Mass [M_sun]':>15} {'Phi (normalized)':>20} {'Plot':>20}")
    print("-" * 60)
    
    for mass, phi in plot_data:
        normalized_phi = phi / max_phi
        bar_length = int(normalized_phi * 30)
        bar = "*" * bar_length
        
        print(f"{mass:>15.2e} {phi:>20.2e} {bar}")
    
    print("-" * 60)

def fit_simple_schechter(mass_centers, mass_function):
    """
    Simple parameter estimation for Schechter function
    (Not a full fit, just rough estimates)
    """
    print("\n" + "="*50)
    print("ROUGH SCHECHTER FUNCTION PARAMETER ESTIMATES")
    print("="*50)
    
    # Filter positive values
    valid_data = [(m, p) for m, p in zip(mass_centers, mass_function) if p > 0]
    
    if not valid_data:
        print("No valid data for fitting")
        return
    
    # Find peak (M* estimate)
    max_phi = max(p for m, p in valid_data)
    M_star_estimate = next(m for m, p in valid_data if p == max_phi)
    
    # phi* estimate
    phi_star_estimate = max_phi
    
    # Simple alpha estimate (slope at low masses)
    # Use first few points to estimate slope
    if len(valid_data) >= 3:
        m1, p1 = valid_data[0]
        m2, p2 = valid_data[1]
        
        if p1 > 0 and p2 > 0 and m1 != m2:
            alpha_estimate = math.log(p2/p1) / math.log(m2/m1)
        else:
            alpha_estimate = -1.0
    else:
        alpha_estimate = -1.0
    
    print(f"φ* ≈ {phi_star_estimate:.2e} Mpc^-3 dex^-1")
    print(f"M* ≈ {M_star_estimate:.2e} M_sun")
    print(f"α  ≈ {alpha_estimate:.2f}")
    print("\nNote: These are rough estimates, not proper fits")

def main():
    """Main analysis function"""
    print("=== Basic Stellar Mass Function Calculator ===\n")
    
    # Load data
    data = load_galaxy_data('Task1_Data_Galaxy_Prop_33.dat')
    
    if not data:
        print("No data available")
        return
    
    # Extract stellar masses (assume column 1, 0-indexed)
    stellar_mass_column = 1 if len(data[0]) >= 2 else 0
    
    # Get stellar masses
    log_stellar_masses = [row[stellar_mass_column] for row in data]
    
    # Check if masses are in log units (typical if values < 15)
    if all(m < 15 for m in log_stellar_masses):
        print("Converting from log10 to linear stellar masses")
        stellar_masses = [10**m for m in log_stellar_masses]
    else:
        stellar_masses = log_stellar_masses
    
    print(f"\nStellar mass statistics:")
    print(f"- Number of galaxies: {len(stellar_masses)}")
    print(f"- Mass range: {min(stellar_masses):.2e} - {max(stellar_masses):.2e} M_sun")
    print(f"- Median mass: {sorted(stellar_masses)[len(stellar_masses)//2]:.2e} M_sun")
    
    # Calculate mass function
    mass_centers, phi, phi_err = calculate_stellar_mass_function(
        stellar_masses, 
        mass_min=1e8, 
        mass_max=1e12,
        n_bins=15,
        survey_volume=1e6  # Mpc^3
    )
    
    # Print results
    print_results(mass_centers, phi, phi_err)
    
    # Create simple plot
    create_simple_plot_data(mass_centers, phi)
    
    # Simple Schechter function estimates
    fit_simple_schechter(mass_centers, phi)
    
    print("\n=== Analysis Complete ===")
    print("For proper plots, install matplotlib and run the full version.")
    
    return mass_centers, phi, phi_err

if __name__ == "__main__":
    mass_centers, phi, phi_err = main()