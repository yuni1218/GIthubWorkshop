# Stellar Mass Function Calculator

This repository contains Python scripts to calculate and plot the stellar mass function from galaxy properties data.

## Overview

The stellar mass function describes the number density of galaxies as a function of their stellar mass. It's a fundamental tool in galaxy evolution studies and provides insights into how galaxies form and evolve over cosmic time.

## Files

- `stellar_mass_function.py` - Comprehensive class-based implementation with advanced features (requires numpy, matplotlib, pandas, scipy, astropy)
- `simple_smf_calculator.py` - Simple, standalone script for basic calculations (requires numpy, matplotlib, pandas, scipy)
- `basic_smf_calculator.py` - Basic version using only standard Python libraries (no external dependencies)
- `requirements.txt` - Required Python packages for full versions
- `example_Task1_Data_Galaxy_Prop_33.dat` - Example data file format
- `README.md` - This file

## Installation

1. Install required packages:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Calculator (No external dependencies)

```bash
python3 basic_smf_calculator.py
```

This script:
- Uses only standard Python libraries
- Provides ASCII-based plotting
- Calculates stellar mass function with error bars
- Gives rough Schechter function parameter estimates
- Perfect for testing and understanding the algorithm

### Simple Calculator (Recommended)

```bash
python simple_smf_calculator.py
```

This script will:
1. Try to load `Task1_Data_Galaxy_Prop_33.dat`
2. If the file is not found, generate realistic sample data
3. Calculate the stellar mass function
4. Create proper plots with error bars
5. Fit a Schechter function to the data

### Advanced Calculator

```bash
python stellar_mass_function.py
```

This provides more advanced features including cosmological calculations and comparison plots.

## Data Format

The scripts expect data in `Task1_Data_Galaxy_Prop_33.dat` with:
- Space-separated values
- Column 1 (index 0): Galaxy ID
- Column 2 (index 1): Stellar mass (log10 or linear units)
- Additional columns: Other properties (redshift, etc.)

## Theory

The stellar mass function is calculated as:

```
Φ(M) = N(M) / (V_survey × Δlog M)
```

Where:
- N(M) = number of galaxies in mass bin M
- V_survey = survey volume in Mpc³
- Δlog M = logarithmic width of mass bin in dex

The Schechter function fit has the form:
```
Φ(M) = φ* × (M/M*)^α × exp(-M/M*)
```

Where:
- φ* = characteristic number density
- M* = characteristic mass
- α = low-mass slope

## Output

The scripts generate:
- `stellar_mass_function.png` - Basic mass function plot
- `stellar_mass_function_with_fit.png` - Plot with Schechter function fit
- Console output with fit parameters and statistics

## Customization

You can modify parameters in the scripts:
- `mass_range`: Stellar mass range for calculation
- `n_bins`: Number of mass bins
- `survey_volume`: Survey volume in Mpc³ (affects normalization)

## Example Output

The code will print results like:
```
Schechter Function Fit:
φ* = 3.45e-03 ± 2.1e-04 Mpc^-3 dex^-1
M* = 2.51e+10 ± 1.2e+09 M_sun
α  = -1.25 ± 0.08
```

## Notes

- If `Task1_Data_Galaxy_Prop_33.dat` is not found, the scripts will generate realistic sample data
- The code automatically detects if stellar masses are in log10 or linear units
- Error bars are calculated using Poisson statistics
- All plots are saved as high-resolution PNG files
