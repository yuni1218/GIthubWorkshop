# Stellar Mass Function Calculation - Complete Solution

## What is the Stellar Mass Function?

The stellar mass function (SMF) is a fundamental tool in galaxy evolution studies that describes the number density of galaxies as a function of their stellar mass. It tells us how many galaxies exist at different mass scales in a given volume of the universe.

### Mathematical Definition

The stellar mass function Φ(M) is defined as:

```
Φ(M) = N(M) / (V_survey × Δlog M)
```

Where:
- **N(M)** = number of galaxies in stellar mass bin M
- **V_survey** = survey volume in Mpc³
- **Δlog M** = logarithmic width of mass bin (in dex)

Units: **[Mpc⁻³ dex⁻¹]** (galaxies per cubic megaparsec per decade in mass)

### Schechter Function

The stellar mass function is often fit with a Schechter function:

```
Φ(M) = φ* × (M/M*)^α × exp(-M/M*)
```

Parameters:
- **φ*** = characteristic number density [Mpc⁻³ dex⁻¹]
- **M*** = characteristic stellar mass [M☉]
- **α** = low-mass slope (typically negative)

## Implementation Details

### 1. Data Loading
The code expects data in `Task1_Data_Galaxy_Prop_33.dat` with columns:
- Column 1: Galaxy ID
- Column 2: Stellar mass (log₁₀ or linear units)
- Additional columns: Other properties (redshift, coordinates, etc.)

### 2. Mass Function Calculation Algorithm

1. **Create logarithmic mass bins**: Evenly spaced in log space from 10⁸ to 10¹² M☉
2. **Count galaxies**: Histogram galaxies into mass bins
3. **Calculate bin centers**: Geometric mean of bin edges
4. **Compute number density**: N(M) / (V_survey × Δlog M)
5. **Calculate errors**: Poisson statistics (√N / (V_survey × Δlog M))

### 3. Key Features

- **Automatic unit detection**: Converts log₁₀ masses to linear if needed
- **Flexible data format**: Handles various column separators and comments
- **Error calculation**: Proper Poisson error propagation
- **Schechter fitting**: Non-linear least squares fitting
- **Visualization**: Multiple plot types with proper formatting

## Code Versions Provided

### 1. `basic_smf_calculator.py` ⭐ RECOMMENDED FOR UNDERSTANDING
**Dependencies**: None (only standard Python libraries)

**Features**:
- Complete stellar mass function calculation
- ASCII-based plotting
- Rough Schechter parameter estimates
- Educational output with detailed explanations
- Works immediately without external packages

**Best for**: Understanding the algorithm, testing, environments without scientific Python packages

### 2. `simple_smf_calculator.py` ⭐ RECOMMENDED FOR RESEARCH
**Dependencies**: numpy, matplotlib, pandas, scipy

**Features**:
- Professional-quality plots
- Proper Schechter function fitting
- Error bars and statistical analysis
- Publication-ready figures
- Handles larger datasets efficiently

**Best for**: Research work, publication figures, detailed analysis

### 3. `stellar_mass_function.py` ⭐ ADVANCED FEATURES
**Dependencies**: numpy, matplotlib, pandas, scipy, astropy

**Features**:
- Object-oriented design
- Cosmological calculations
- Advanced plotting options
- Literature comparison capabilities
- Comprehensive error handling

**Best for**: Advanced research, cosmological studies, large surveys

## Running the Code

### Quick Start (No Dependencies)
```bash
python3 basic_smf_calculator.py
```

### With Scientific Python
```bash
# Install dependencies
pip install numpy matplotlib pandas scipy

# Run analysis
python simple_smf_calculator.py
```

### With Your Data
1. Place your data file as `Task1_Data_Galaxy_Prop_33.dat`
2. Run any of the scripts
3. Check output plots and printed results

## Expected Output

### Console Output Example
```
=== Stellar Mass Function Calculator ===

Successfully loaded 5000 galaxies from Task1_Data_Galaxy_Prop_33.dat
Converting from log10 to linear stellar masses

Stellar mass statistics:
- Number of galaxies: 5000
- Mass range: 1.32e+08 - 7.94e+11 M_sun
- Median mass: 1.58e+10 M_sun

Calculated stellar mass function:
- Total galaxies in range: 4987
- Mass range: 1.0e+08 - 1.0e+12 M_sun
- Survey volume: 1.0e+06 Mpc^3

Schechter Function Fit:
φ* = 3.45e-03 ± 2.1e-04 Mpc^-3 dex^-1
M* = 2.51e+10 ± 1.2e+09 M_sun
α  = -1.25 ± 0.08
```

### Generated Files
- `stellar_mass_function.png` - Basic mass function plot
- `stellar_mass_function_with_fit.png` - Plot with Schechter fit
- `stellar_mass_function_comparison.png` - Advanced comparison plot

## Scientific Context

### Typical Values
- **φ*** ≈ 10⁻³ to 10⁻² Mpc⁻³ dex⁻¹
- **M*** ≈ 10¹⁰·⁵ to 10¹¹ M☉
- **α** ≈ -1.0 to -1.5

### Physical Interpretation
- **High-mass end**: Exponential cutoff due to feedback processes
- **Low-mass end**: Power-law behavior, slope α indicates efficiency of galaxy formation
- **Knee**: Characteristic mass M* represents transition between regimes

### Applications
- Galaxy formation and evolution studies
- Dark matter halo mass functions
- Cosmological parameter constraints
- Survey completeness testing

## Customization Options

### Parameters to Adjust
```python
# Mass range
mass_range = (1e8, 1e12)  # M_sun

# Number of bins
n_bins = 20

# Survey volume (affects normalization)
survey_volume = 1e6  # Mpc^3

# Mass column (if different from column 1)
stellar_mass_column = 1  # 0-indexed
```

### Adding Features
- **Completeness corrections**: Multiply by detection efficiency
- **Cosmic variance**: Add systematic uncertainties
- **Redshift evolution**: Calculate for redshift bins
- **Environment dependence**: Split by galaxy environment

## Troubleshooting

### Common Issues
1. **"No module named 'numpy'"**: Use `basic_smf_calculator.py` instead
2. **Empty plots**: Check mass range and survey volume
3. **Fitting failures**: Increase number of galaxies or adjust mass range
4. **File not found**: Code will generate sample data automatically

### Data Format Problems
- Ensure space-separated columns
- Remove header lines (or prefix with #)
- Check for missing values or text in numeric columns

## Theory Background

### Why Logarithmic Bins?
Galaxy masses span many orders of magnitude (10⁸ to 10¹² M☉), so logarithmic binning provides equal statistical weight across the mass range.

### Poisson Errors
Galaxy counts follow Poisson statistics, so the error in each bin is √N, where N is the number of galaxies in that bin.

### Survey Volume
The survey volume depends on the geometry and redshift range of your survey. For redshift z, the comoving volume element is proportional to D_C²(z) dz, where D_C is the comoving distance.

### Schechter Function Origin
The Schechter function arises from the Press-Schechter formalism for dark matter halo formation, modified by baryonic physics that affects star formation efficiency.

---

**This implementation provides a complete, research-grade stellar mass function calculator suitable for educational use, research applications, and publication-quality analysis.**