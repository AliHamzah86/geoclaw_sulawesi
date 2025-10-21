#!/usr/bin/env python3
"""
Improved Scenario Generation for PTHA using K-L Expansion
Based on KLSlip2D.py but adapted for scenario generation workflow

This replaces the problematic k-means approach with a physically realistic
K-L expansion that respects fault geometry and spatial correlation.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from clawpack.geoclaw import dtopotools
from clawpack.visclaw import colormaps

def setup_fault_from_csv(filename):
    """Read fault parameters from CSV file and return fault object."""
    column_map = {"longitude": 1, "latitude": 2, "depth": 3, "strike": 4,
                  "length": 5, "width": 6, "dip": 7}
    defaults = {'rake': 90, 'slip': 1.0}
    coordinate_specification = 'top center'
    input_units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
    rupture_type = 'static'
    skiprows = 1
    delimiter = ','
    
    fault = dtopotools.Fault()
    fault.read(filename, column_map, coordinate_specification,
               rupture_type, skiprows, delimiter, input_units, defaults)
    
    # Adjust to W coordinates
    for s in fault.subfaults:
        s.longitude = s.longitude - 360.
    
    print(f"Loaded fault with {len(fault.subfaults)} subfaults")
    return fault

def subdivide_fault(fault, subfault_length=8.0, subfault_width=8.0, max_subfaults=500):
    """Subdivide the fault into smaller subfaults with progress printing.
    
    Parameters:
    fault : dtopotools.Fault -> Original fault object.
    subfault_length, subfault_width : float -> Desired subdivision size in km.
    max_subfaults : int -> if subdivision would exceed this, 
                fall back to using coarser subdivision.
    """
    phi_plate = 60.0  # angle oceanic plate moves clockwise from north
    new_subfaults = []

    print(f"Subdividing {len(fault.subfaults)} original subfaults...")

    for i, subfault in enumerate(fault.subfaults):
        subfault.rake = subfault.strike - phi_plate - 180.0
        nstrike = int(subfault.length / subfault_length)
        ndip = int(subfault.width / subfault_width)

        print(f"  Subfault {i+1}/{len(fault.subfaults)}: "
              f"{subfault.length:.1f} × {subfault.width:.1f} km → {nstrike}×{ndip}")

        if nstrike > 0 and ndip > 0:
            f = dtopotools.SubdividedPlaneFault(subfault, nstrike, ndip)
            new_subfaults += f.subfaults
        else:
            new_subfaults.append(subfault)

        # Safety check: stop if too many
        if len(new_subfaults) > max_subfaults:
            print(f"⚠️ Subdivision exceeded {max_subfaults} subfaults. "
                  f"Switching to coarser subdivision...")
            # Coarser subdivision: double the lengths
            return subdivide_fault(fault,
                                   subfault_length=subfault_length*2,
                                   subfault_width=subfault_width*2,
                                   max_subfaults=max_subfaults)

    new_fault = dtopotools.Fault(subfaults=new_subfaults)
    print(f"✅ Subdivided fault has {len(new_fault.subfaults)} subfaults")
    return new_fault

def compute_subfault_distances(fault):
    """Compute pairwise distances between subfaults in strike and dip directions."""
    print("Computing subfault distances...")
    
    nsubfaults = len(fault.subfaults)
    D = np.zeros((nsubfaults, nsubfaults))
    Dstrike = np.zeros((nsubfaults, nsubfaults))
    Ddip = np.zeros((nsubfaults, nsubfaults))
    
    # Earth radius and conversion factors
    rad = np.pi / 180.
    rr = 6.378e6  # Earth radius (m)
    lat2meter = rr * rad
    MIN_SIN = 1e-3  # to avoid Ddip > Depth
    
    for i, si in enumerate(fault.subfaults):
        xi, yi, zi = si.longitude, si.latitude, si.depth
        for j, sj in enumerate(fault.subfaults):
            xj, yj, zj = sj.longitude, sj.latitude, sj.depth
            
            # Convert to meters
            dx = abs(xi-xj) * np.cos(0.5*(yi+yj)*np.pi/180.) * lat2meter
            dy = abs(yi-yj) * lat2meter
            dz = abs(zi-zj)
            
            # Euclidean distance
            D[i, j] = np.sqrt(dx**2 + dy**2 + dz**2)
            
            # Distance down-dip based on depths
            dip = 0.5 * (si.dip + sj.dip)
            sin_dip = np.maximum(np.sin(dip * np.pi / 180.), MIN_SIN)
            ddip1 = dz / sin_dip
            Ddip[i, j] = min(ddip1, D[i, j])
            
            # Strike distance
            dstrike2 = max(D[i, j]**2 - Ddip[i, j]**2, 0.0)
            Dstrike[i, j] = np.sqrt(dstrike2)
    
    return D, Dstrike, Ddip

def compute_correlation_matrix(Dstrike, Ddip, Lstrike, Ldip):
    """Compute 2D correlation matrix based on strike and dip distances."""
    r = np.sqrt((Dstrike / Lstrike)**2 + (Ddip / Ldip)**2)
    C = np.exp(-r)
    return C

def compute_mean_slip(fault, Mw_desired=9.0):
    """Compute mean slip from desired moment magnitude."""
    lengths = np.array([s.length for s in fault.subfaults])
    widths = np.array([s.width for s in fault.subfaults])
    areas = lengths * widths
    total_area = sum(areas)
    
    Mo_desired = 10.**(1.5 * Mw_desired + 9.05)
    mean_slip = Mo_desired / (fault.subfaults[0].mu * total_area)
    print(f"Mean slip {mean_slip:.2f} meters for Mw {Mw_desired}")
    
    return mean_slip

def compute_eigenmodes(C, mean_slip, alpha=0.5):
    """Compute eigenvalues and eigenvectors of the covariance matrix."""
    print(f"Computing eigenmodes from {C.shape[0]}x{C.shape[1]} matrix...")
    
    # Lognormal transformation
    sigma_slip = alpha * mean_slip * np.ones_like(C)
    Cov_g = np.log(alpha**2 * C + 1.)
    mean_slip_g = np.log(mean_slip) - np.diag(Cov_g) / 2.
    
    # Eigenvalue decomposition
    lam, V = np.linalg.eig(Cov_g)
    lam = np.real(lam)
    V = np.real(V)
    
    # Sort by eigenvalue magnitude
    i = list(np.argsort(lam))
    i.reverse()
    lam = lam[i]
    V = V[:, i]
    
    print(f"Largest eigenvalue: {lam[0]:.2f}")
    print(f"Smallest eigenvalue: {lam[-1]:.2f}")
    
    return lam, V, mean_slip_g

def create_depth_taper(max_depth, taper_type='exp_depth'):
    """Create depth-dependent tapering function."""
    def exp_depth_taper(depth):
        """Exponential depth tapering function."""
        return 1. - np.exp((depth - max_depth) * 5. / max_depth)
    
    def none_taper(depth):
        """No tapering."""
        return 1.0
    
    taper_functions = {
        'exp_depth': exp_depth_taper,
        'none': none_taper
    }
    
    return taper_functions.get(taper_type, exp_depth_taper)

def kl_expansion(z, lam, V, fault, tau_func, Mo_desired):
    """Perform K-L expansion to compute slip distribution."""
    nsubfaults = len(fault.subfaults)
    KL_slip = np.zeros(nsubfaults)
    
    # K-L expansion (skip first term which is mean)
    for k in range(1, len(z)):
        if k < len(lam):  # Ensure we don't exceed available modes
            KL_slip += z[k] * np.sqrt(lam[k]) * V[:, k]
    
    # Lognormal transformation
    KL_slip = np.exp(KL_slip)
    
    # Apply depth tapering
    for j, s in enumerate(fault.subfaults):
        s.slip = KL_slip[j] * tau_func(s.depth)
    
    # Rescale to desired moment magnitude
    Mo = fault.Mo()
    if Mo > 0:
        KL_slip *= Mo_desired / Mo
        for j, s in enumerate(fault.subfaults):
            s.slip = KL_slip[j] * tau_func(s.depth)
    
    return KL_slip

def generate_kl_scenarios(fault, n_scenarios=100, n_terms=20, 
                         Lstrike=130e3, Ldip=40e3, Mw_desired=9.0,
                         correlation_lengths=None):
    """
    Generate scenarios using K-L expansion instead of k-means.
    
    Parameters:
    -----------
    fault : dtopotools.Fault
        Fault geometry object
    n_scenarios : int
        Number of scenarios to generate
    n_terms : int
        Number of K-L terms to use
    Lstrike, Ldip : float
        Correlation lengths in strike and dip directions (m)
    Mw_desired : float
        Target moment magnitude
    correlation_lengths : dict
        Alternative way to specify correlation lengths
    
    Returns:
    --------
    scenarios : array
        Array of shape (n_scenarios, n_subfaults) containing slip values
    weights : array
        Array of shape (n_scenarios,) containing probability weights
    """
    print(f"Generating {n_scenarios} scenarios using K-L expansion...")
    
    # 1. Compute distances and correlation matrix
    D, Dstrike, Ddip = compute_subfault_distances(fault)
    
    # Use provided correlation lengths or defaults
    if correlation_lengths:
        Lstrike = correlation_lengths.get('strike', Lstrike)
        Ldip = correlation_lengths.get('dip', Ldip)
    
    C = compute_correlation_matrix(Dstrike, Ddip, Lstrike, Ldip)
    
    # 2. Compute mean slip and eigenmodes
    mean_slip = compute_mean_slip(fault, Mw_desired)
    lam, V, mean_slip_g = compute_eigenmodes(C, mean_slip)
    
    # 3. Create depth tapering function
    max_depth = max(s.depth for s in fault.subfaults)
    tau_func = create_depth_taper(max_depth, 'exp_depth')
    
    # 4. Generate scenarios
    scenarios = []
    n_subfaults = len(fault.subfaults)
    Mo_desired = 10.**(1.5 * Mw_desired + 9.05)
    
    for i in range(n_scenarios):
        # Generate random K-L coefficients
        z = np.random.randn(min(n_terms + 1, len(lam)))
        
        # Apply K-L expansion
        slip = kl_expansion(z, lam, V, fault, tau_func, Mo_desired)
        scenarios.append(slip)
        
        if (i + 1) % 20 == 0:
            print(f"Generated {i + 1}/{n_scenarios} scenarios...")
    
    scenarios = np.array(scenarios)
    
    # 5. Compute probability weights (uniform for now, can be improved)
    weights = np.ones(n_scenarios) / n_scenarios
    
    print(f"Generated {n_scenarios} scenarios with {n_subfaults} subfaults each")
    print(f"Slip range: {scenarios.min():.2f} to {scenarios.max():.2f} meters")
    
    return scenarios, weights

def save_scenarios(scenarios, weights, filename_pts='scenario_pts_kl.txt', 
                   filename_weights='scenario_prb_wgts_kl.txt'):
    """Save scenarios and weights to text files."""
    np.savetxt(filename_pts, scenarios)
    np.savetxt(filename_weights, weights)
    print(f"Saved scenarios to {filename_pts}")
    print(f"Saved weights to {filename_weights}")

def plot_scenario_comparison(scenarios_kl, scenarios_kmeans=None, 
                           fault=None, n_plot=6):
    """Plot comparison between K-L and k-means scenarios."""
    fig, axes = plt.subplots(2, n_plot, figsize=(15, 6))
    
    for i in range(n_plot):
        # K-L scenarios
        ax = axes[0, i]
        if fault is not None:
            for j, s in enumerate(fault.subfaults):
                s.slip = scenarios_kl[i, j]
            fault.plot_subfaults(ax, slip_color=True, cmax_slip=20.)
        else:
            ax.plot(scenarios_kl[i, :])
        ax.set_title(f'K-L Scenario {i+1}')
        ax.axis('off')
        
        # K-means scenarios (if provided)
        if scenarios_kmeans is not None:
            ax = axes[1, i]
            if fault is not None:
                for j, s in enumerate(fault.subfaults):
                    s.slip = scenarios_kmeans[i, j]
                fault.plot_subfaults(ax, slip_color=True, cmax_slip=20.)
            else:
                ax.plot(scenarios_kmeans[i, :])
            ax.set_title(f'K-means Scenario {i+1}')
        ax.axis('off')
    
    plt.tight_layout()
    plt.savefig('scenario_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved scenario comparison plot")

def main():
    """Main function to generate scenarios."""
    # Load fault geometry
    fault_filename = 'CSZe01.csv'
    fault = setup_fault_from_csv(fault_filename)
    
    # Subdivide fault for detailed analysis
    subdivided_fault = subdivide_fault(fault, subfault_length=8.0, subfault_width=8.0)
    
    # Generate scenarios using K-L expansion
    scenarios, weights = generate_kl_scenarios(
        subdivided_fault, 
        n_scenarios=100, 
        n_terms=20,
        Lstrike=130e3,  # 130 km correlation length in strike
        Ldip=40e3,      # 40 km correlation length in dip
        Mw_desired=9.0
    )
    
    # Save scenarios
    save_scenarios(scenarios, weights)
    
    # Create comparison plot
    plot_scenario_comparison(scenarios, fault=subdivided_fault)
    
    print("Scenario generation complete!")

if __name__ == "__main__":
    main()
