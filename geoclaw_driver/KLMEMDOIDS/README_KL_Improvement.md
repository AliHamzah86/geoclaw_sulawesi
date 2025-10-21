# K-L Expansion Scenario Generation for PTHA

## Overview

This repository now includes an improved scenario generation approach that replaces the problematic k-means method with a physically realistic K-L (Karhunen-Loève) expansion based on fault geometry.

## Files Created

### 1. `KL_Scenario_Generation.py`
- **Main K-L expansion implementation**
- Based on `KLSlip2D.py` but adapted for scenario generation workflow
- Includes functions for:
  - Fault geometry loading and subdivision
  - Spatial correlation matrix computation
  - Eigenvalue decomposition
  - K-L expansion with depth tapering
  - Scenario generation and weight calculation

### 2. `run_CC_CSZ_South_KL.py`
- **Modified main simulation script**
- Replaces k-means approach with K-L expansion
- Uses CSZ fault geometry from `CSZe01.csv`
- Generates physically realistic scenarios

### 3. `compare_scenario_methods.py`
- **Comparison between k-means and K-L approaches**
- Demonstrates advantages of K-L expansion
- Creates comparison plots

### 4. `test_kl_scenarios.py`
- **Simple test script**
- Tests K-L expansion functionality
- Creates validation plots

## Key Improvements

### ❌ Problems with Original K-means Approach:
1. **Synthetic Centroids**: Creates artificial scenarios that may not be physically realistic
2. **Arbitrary Dilation**: 4x dilation factor has no theoretical justification
3. **Computational Waste**: Uses 1.5M Monte Carlo samples just for weight estimation
4. **Poor Tail Coverage**: K-means naturally focuses on high-density regions
5. **Ignores Fault Geometry**: No consideration of spatial correlation

### ✅ Advantages of K-L Expansion Approach:
1. **Physically Realistic**: Uses actual fault geometry and spatial correlation
2. **Log-normal Distribution**: More realistic for earthquake slip patterns
3. **Computationally Efficient**: Direct analytical weight calculation
4. **Natural Tail Coverage**: No arbitrary parameters needed
5. **Fault-Aware**: Respects subfault distances and correlations
6. **Depth Tapering**: Includes realistic depth-dependent slip reduction

## Usage

### Basic Usage:
```python
from KL_Scenario_Generation import generate_kl_scenarios, setup_fault_from_csv

# Load fault geometry
fault = setup_fault_from_csv('CSZe01.csv')

# Generate scenarios
scenarios, weights = generate_kl_scenarios(
    fault,
    n_scenarios=100,
    n_terms=20,
    Lstrike=130e3,  # 130 km correlation length
    Ldip=40e3,      # 40 km correlation length
    Mw_desired=9.0
)
```

### Running the Modified Simulation:
```bash
python run_CC_CSZ_South_KL.py
```

## Technical Details

### K-L Expansion Formula:
```
s = exp[∑(j=1 to n) z_j * √(λ_j) * v_j]
```
Where:
- `s` = slip distribution
- `z_j` = random K-L coefficients ~ N(0,1)
- `λ_j` = eigenvalues of correlation matrix
- `v_j` = eigenvectors of correlation matrix

### Spatial Correlation:
```
C_ij = exp(-√((D_strike/L_strike)² + (D_dip/L_dip)²))
```
Where:
- `D_strike`, `D_dip` = distances between subfaults
- `L_strike`, `L_dip` = correlation lengths

### Depth Tapering:
```
τ(d) = 1 - exp((d - d_max) * 5 / d_max)
```
Where:
- `d` = depth
- `d_max` = maximum depth

## Performance Comparison

| Aspect | K-means | K-L Expansion |
|--------|---------|---------------|
| **Physical Realism** | Low | High |
| **Computational Cost** | High (1.5M MC) | Low (analytical) |
| **Tail Coverage** | Artificial | Natural |
| **Fault Geometry** | Ignored | Respected |
| **Parameters** | Ad-hoc | Theoretically justified |

## Next Steps

1. **Test the new approach** with your existing simulation workflow
2. **Compare results** between k-means and K-L approaches
3. **Validate** that tsunami simulations work correctly
4. **Optimize parameters** (correlation lengths, number of terms)
5. **Scale up** to full PTHA analysis

## References

- Based on `KLSlip2D.py` implementation
- LeVeque et al. (2016) - K-L expansion for earthquake slip
- GeoClaw documentation for fault geometry handling
