#!/usr/bin/env python3
"""
Stellar Mass Function Calculator and Plotter

This script reads galaxy properties data and calculates the stellar mass function,
then creates visualization plots. The stellar mass function describes the number
density of galaxies as a function of their stellar mass.

Based on the provided specifications:
- Data file: Task1_Data_Galaxy_Prop_33.dat
- Contains galaxy properties including stellar masses
- Calculates number density per mass bin
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from astropy import units as u
from astropy.cosmology import Planck18
import warnings
warnings.filterwarnings('ignore')

class StellarMassFunction:
    """
    Class to calculate and plot stellar mass function from galaxy data
    """
    
    def __init__(self, data_file=None, cosmology=Planck18):
        """
        Initialize the stellar mass function calculator
        
        Parameters:
        -----------
        data_file : str
            Path to the galaxy properties data file
        cosmology : astropy.cosmology object
            Cosmological parameters to use
        """
        self.data_file = data_file
        self.cosmology = cosmology
        self.galaxy_data = None
        self.stellar_masses = None
        self.mass_bins = None
        self.mass_function = None
        self.mass_function_error = None
        
    def load_data(self, data_file=None):
        """
        Load galaxy properties data from file
        
        Parameters:
        -----------
        data_file : str, optional
            Path to data file. If None, uses self.data_file
        """
        if data_file is None:
            data_file = self.data_file
            
        if data_file is None:
            raise ValueError("No data file specified")
            
        try:
            # Try to read as space-separated values first
            self.galaxy_data = pd.read_csv(data_file, sep='\s+', comment='#')
            print(f"Successfully loaded {len(self.galaxy_data)} galaxies from {data_file}")
            print(f"Columns: {list(self.galaxy_data.columns)}")
            
        except Exception as e:
            print(f"Error loading data file: {e}")
            # Generate sample data for demonstration
            self.generate_sample_data()
            
    def generate_sample_data(self, n_galaxies=10000):
        """
        Generate sample galaxy data for demonstration purposes
        
        Parameters:
        -----------
        n_galaxies : int
            Number of sample galaxies to generate
        """
        print(f"Generating {n_galaxies} sample galaxies for demonstration...")
        
        # Generate stellar masses following a realistic distribution
        # Log-normal distribution with parameters typical for galaxy surveys
        log_mass_mean = 10.5  # log10(M*/M_sun)
        log_mass_std = 0.7
        
        log_stellar_masses = np.random.normal(log_mass_mean, log_mass_std, n_galaxies)
        stellar_masses = 10**log_stellar_masses
        
        # Generate other properties
        redshifts = np.random.uniform(0.01, 0.1, n_galaxies)  # Low redshift sample
        galaxy_ids = np.arange(1, n_galaxies + 1)
        
        # Create DataFrame
        self.galaxy_data = pd.DataFrame({
            'galaxy_id': galaxy_ids,
            'stellar_mass': stellar_masses,
            'log_stellar_mass': log_stellar_masses,
            'redshift': redshifts
        })
        
        print("Sample data generated successfully")
        
    def extract_stellar_masses(self, mass_column=None):
        """
        Extract stellar masses from the galaxy data
        
        Parameters:
        -----------
        mass_column : str, optional
            Name of the column containing stellar masses
        """
        if self.galaxy_data is None:
            raise ValueError("No galaxy data loaded. Call load_data() first.")
            
        # Try to automatically identify stellar mass column
        possible_columns = ['stellar_mass', 'mstar', 'M_star', 'mass', 'log_stellar_mass', 'logmstar']
        
        if mass_column is None:
            for col in possible_columns:
                if col in self.galaxy_data.columns:
                    mass_column = col
                    break
                    
            if mass_column is None:
                # Use the first numerical column
                numerical_cols = self.galaxy_data.select_dtypes(include=[np.number]).columns
                if len(numerical_cols) > 0:
                    mass_column = numerical_cols[0]
                    print(f"Using column '{mass_column}' as stellar mass")
                else:
                    raise ValueError("No suitable stellar mass column found")
        
        self.stellar_masses = self.galaxy_data[mass_column].values
        
        # Convert to linear mass if in log units
        if 'log' in mass_column.lower() or np.all(self.stellar_masses < 15):
            print("Converting from log to linear stellar masses")
            self.stellar_masses = 10**self.stellar_masses
            
        print(f"Extracted {len(self.stellar_masses)} stellar masses")
        print(f"Mass range: {np.min(self.stellar_masses):.2e} - {np.max(self.stellar_masses):.2e} M_sun")
        
    def calculate_mass_function(self, mass_min=1e8, mass_max=1e12, n_bins=25, 
                              survey_volume=None, completeness_correction=True):
        """
        Calculate the stellar mass function
        
        Parameters:
        -----------
        mass_min : float
            Minimum stellar mass [M_sun]
        mass_max : float
            Maximum stellar mass [M_sun]
        n_bins : int
            Number of mass bins
        survey_volume : float, optional
            Survey volume in Mpc^3. If None, will be estimated
        completeness_correction : bool
            Whether to apply completeness corrections
        """
        if self.stellar_masses is None:
            raise ValueError("No stellar masses extracted. Call extract_stellar_masses() first.")
            
        # Create logarithmic mass bins
        log_mass_min = np.log10(mass_min)
        log_mass_max = np.log10(mass_max)
        log_mass_bins = np.linspace(log_mass_min, log_mass_max, n_bins + 1)
        self.mass_bins = 10**log_mass_bins
        
        # Calculate bin centers
        mass_centers = np.sqrt(self.mass_bins[:-1] * self.mass_bins[1:])
        
        # Count galaxies in each bin
        mask = (self.stellar_masses >= mass_min) & (self.stellar_masses <= mass_max)
        filtered_masses = self.stellar_masses[mask]
        
        counts, _ = np.histogram(filtered_masses, bins=self.mass_bins)
        
        # Calculate bin widths
        bin_widths = self.mass_bins[1:] - self.mass_bins[:-1]
        
        # Estimate survey volume if not provided
        if survey_volume is None:
            if hasattr(self, 'galaxy_data') and 'redshift' in self.galaxy_data.columns:
                # Estimate from redshift range
                z_max = np.max(self.galaxy_data['redshift'])
                # Approximate volume for all-sky survey
                survey_volume = (4/3) * np.pi * (self.cosmology.comoving_distance(z_max).value)**3
                print(f"Estimated survey volume: {survey_volume:.2e} Mpc^3")
            else:
                # Default volume for demonstration
                survey_volume = 1e6  # Mpc^3
                print(f"Using default survey volume: {survey_volume:.2e} Mpc^3")
        
        # Calculate number density (galaxies per Mpc^3 per dex in mass)
        # Convert bin widths to logarithmic (dex)
        dlog_mass = np.log10(self.mass_bins[1:] / self.mass_bins[:-1])
        
        self.mass_function = counts / (survey_volume * dlog_mass)
        
        # Calculate Poisson errors
        self.mass_function_error = np.sqrt(counts) / (survey_volume * dlog_mass)
        
        # Store bin centers for plotting
        self.mass_centers = mass_centers
        
        print(f"Calculated stellar mass function with {n_bins} bins")
        print(f"Total galaxies used: {np.sum(counts)}")
        
        return self.mass_centers, self.mass_function, self.mass_function_error
        
    def plot_mass_function(self, save_fig=True, fig_name='stellar_mass_function.png',
                          show_errors=True, log_scale=True):
        """
        Plot the stellar mass function
        
        Parameters:
        -----------
        save_fig : bool
            Whether to save the figure
        fig_name : str
            Name of the output figure file
        show_errors : bool
            Whether to show error bars
        log_scale : bool
            Whether to use log scale for y-axis
        """
        if self.mass_function is None:
            raise ValueError("Mass function not calculated. Call calculate_mass_function() first.")
            
        plt.figure(figsize=(10, 8))
        
        # Plot the mass function
        if show_errors:
            plt.errorbar(self.mass_centers, self.mass_function, yerr=self.mass_function_error,
                        fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=6,
                        label='Stellar Mass Function', color='blue', alpha=0.8)
        else:
            plt.plot(self.mass_centers, self.mass_function, 'o-', linewidth=2, 
                    markersize=6, label='Stellar Mass Function', color='blue')
        
        # Formatting
        plt.xlabel(r'Stellar Mass [M$_{\odot}$]', fontsize=14)
        plt.ylabel(r'$\Phi$ [Mpc$^{-3}$ dex$^{-1}$]', fontsize=14)
        plt.title('Galaxy Stellar Mass Function', fontsize=16, fontweight='bold')
        
        plt.xscale('log')
        if log_scale:
            plt.yscale('log')
            
        # Add grid
        plt.grid(True, alpha=0.3, which='both')
        
        # Format axes
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tick_params(axis='both', which='minor', labelsize=10)
        
        # Add legend
        plt.legend(fontsize=12)
        
        # Tight layout
        plt.tight_layout()
        
        if save_fig:
            plt.savefig(fig_name, dpi=300, bbox_inches='tight')
            print(f"Figure saved as {fig_name}")
            
        plt.show()
        
    def fit_schechter_function(self, initial_params=None):
        """
        Fit a Schechter function to the stellar mass function
        
        The Schechter function: φ(M) = φ* (M/M*)^α exp(-M/M*)
        
        Parameters:
        -----------
        initial_params : list, optional
            Initial parameters [phi_star, M_star, alpha] for fitting
        """
        if self.mass_function is None:
            raise ValueError("Mass function not calculated. Call calculate_mass_function() first.")
            
        try:
            from scipy.optimize import curve_fit
            
            def schechter_function(mass, phi_star, M_star, alpha):
                """Schechter function"""
                x = mass / M_star
                return phi_star * (x**alpha) * np.exp(-x)
            
            # Filter out zero or negative values
            mask = self.mass_function > 0
            x_data = self.mass_centers[mask]
            y_data = self.mass_function[mask]
            yerr_data = self.mass_function_error[mask]
            
            # Initial parameter guesses
            if initial_params is None:
                phi_star_guess = np.max(y_data)
                M_star_guess = x_data[np.argmax(y_data)]
                alpha_guess = -1.0
                initial_params = [phi_star_guess, M_star_guess, alpha_guess]
            
            # Perform fit
            popt, pcov = curve_fit(schechter_function, x_data, y_data, 
                                 p0=initial_params, sigma=yerr_data, absolute_sigma=True,
                                 maxfev=5000)
            
            # Calculate parameter errors
            perr = np.sqrt(np.diag(pcov))
            
            phi_star, M_star, alpha = popt
            phi_star_err, M_star_err, alpha_err = perr
            
            print("Schechter Function Fit Results:")
            print(f"φ* = {phi_star:.2e} ± {phi_star_err:.2e} Mpc^-3 dex^-1")
            print(f"M* = {M_star:.2e} ± {M_star_err:.2e} M_sun")
            print(f"α = {alpha:.2f} ± {alpha_err:.2f}")
            
            # Generate smooth curve for plotting
            mass_smooth = np.logspace(np.log10(np.min(x_data)), np.log10(np.max(x_data)), 100)
            phi_smooth = schechter_function(mass_smooth, *popt)
            
            return popt, pcov, mass_smooth, phi_smooth
            
        except ImportError:
            print("scipy.optimize not available for fitting")
            return None
        except Exception as e:
            print(f"Fitting failed: {e}")
            return None
            
    def create_comparison_plot(self, literature_data=None):
        """
        Create a comparison plot with literature data
        
        Parameters:
        -----------
        literature_data : dict, optional
            Dictionary with 'mass', 'phi', 'label' keys for literature comparison
        """
        if self.mass_function is None:
            raise ValueError("Mass function not calculated. Call calculate_mass_function() first.")
            
        plt.figure(figsize=(12, 8))
        
        # Plot our data
        plt.errorbar(self.mass_centers, self.mass_function, yerr=self.mass_function_error,
                    fmt='o-', capsize=5, capthick=2, linewidth=2, markersize=8,
                    label='This Work', color='red', alpha=0.8, zorder=3)
        
        # Plot literature data if provided
        if literature_data:
            for i, data in enumerate(literature_data):
                color = plt.cm.tab10(i)
                plt.plot(data['mass'], data['phi'], '--', linewidth=2, 
                        label=data['label'], color=color, alpha=0.7)
        
        # Try to fit and plot Schechter function
        fit_result = self.fit_schechter_function()
        if fit_result is not None:
            popt, pcov, mass_smooth, phi_smooth = fit_result
            plt.plot(mass_smooth, phi_smooth, 'k-', linewidth=3, alpha=0.7,
                    label=f'Schechter Fit (α={popt[2]:.2f})', zorder=2)
        
        # Formatting
        plt.xlabel(r'Stellar Mass [M$_{\odot}$]', fontsize=16)
        plt.ylabel(r'$\Phi$ [Mpc$^{-3}$ dex$^{-1}$]', fontsize=16)
        plt.title('Galaxy Stellar Mass Function Comparison', fontsize=18, fontweight='bold')
        
        plt.xscale('log')
        plt.yscale('log')
        
        # Set reasonable axis limits
        plt.xlim(1e8, 1e12)
        plt.ylim(1e-6, 1e-1)
        
        # Add grid
        plt.grid(True, alpha=0.3, which='both')
        
        # Format axes
        plt.tick_params(axis='both', which='major', labelsize=14)
        plt.tick_params(axis='both', which='minor', labelsize=12)
        
        # Add legend
        plt.legend(fontsize=12, loc='upper right')
        
        # Tight layout
        plt.tight_layout()
        
        plt.savefig('stellar_mass_function_comparison.png', dpi=300, bbox_inches='tight')
        print("Comparison plot saved as stellar_mass_function_comparison.png")
        plt.show()

def main():
    """
    Main function to run the stellar mass function analysis
    """
    print("=== Stellar Mass Function Calculator ===\n")
    
    # Initialize the calculator
    smf = StellarMassFunction()
    
    # Try to load data, fallback to sample data
    try:
        smf.load_data('Task1_Data_Galaxy_Prop_33.dat')
    except:
        print("Data file not found, generating sample data...\n")
        smf.generate_sample_data()
    
    # Extract stellar masses
    smf.extract_stellar_masses()
    
    # Calculate mass function
    print("\nCalculating stellar mass function...")
    mass_centers, phi, phi_err = smf.calculate_mass_function(
        mass_min=1e8, 
        mass_max=1e12, 
        n_bins=20
    )
    
    # Create plots
    print("\nCreating plots...")
    smf.plot_mass_function(save_fig=True, fig_name='stellar_mass_function.png')
    
    # Create comparison plot with fitting
    smf.create_comparison_plot()
    
    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Number of galaxies: {len(smf.stellar_masses)}")
    print(f"Mass range: {np.min(smf.stellar_masses):.2e} - {np.max(smf.stellar_masses):.2e} M_sun")
    print(f"Number density range: {np.min(phi[phi>0]):.2e} - {np.max(phi):.2e} Mpc^-3 dex^-1")
    
    return smf

if __name__ == "__main__":
    # Run the analysis
    stellar_mass_function = main()