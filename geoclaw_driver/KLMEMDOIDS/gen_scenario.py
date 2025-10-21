
from clawpack.geoclaw import dtopotools
import numpy as np
import matplotlib.pyplot as plt


def compute_correlation_matrix(Dstrike, Ddip, Lstrike, Ldip):
    """Compute 2D correlation matrix."""
    r = np.sqrt((Dstrike / Lstrike)**2 + (Ddip / Ldip)**2)
    C = np.exp(-r)
    return C


def compute_eigenmodes(C, mean_slip):
    """Compute eigenvalues and eigenvectors of the covariance matrix."""
    # Lognormal transformation
    Cov_g = np.log(alpha**2 * C + 1.)
    lam, V = np.linalg.eig(Cov_g)
    return lam, V, mean_slip_g

def KL_expansion(z, lam, V, new_fault, tau, Mo_desired):
    """Perform 2D K-L expansion to compute slip."""
    for k in range(1, len(z)):
        KL_slip += z[k] * np.sqrt(lam[k]) * V[:, k]
    KL_slip = np.exp(KL_slip)  # Log-normal
    # Apply depth tapering
    for j, s in enumerate(new_fault.subfaults):
        s.slip = KL_slip[j] * tau(s.depth)


def exp_depth(x):
    tau_depth = lambda d: 1. - np.exp((d - max_depth) * 20 / max_depth)
    return tau_depth(depth0 + x * np.sin(dip * np.pi / 180.))



def setup_fault_from_csv(filename: str) -> dtopotools.Fault:
    """
    Read fault parameters from a CSV file (e.g., 'CSZe01.csv') and return a Fault.
    Assumes columns:
      longitude, latitude, depth, strike, length, width, dip
    Units:
      slip: m, depth/length/width: km
    Coordinates specified at subfault top center.
    """
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

    # Convert longitudes from [0,360) to [-360,0) if needed
    for s in fault.subfaults:
        s.longitude = s.longitude - 360.0

    return fault


def subdivide_fault(fault: dtopotools.Fault,
                    subfault_length: float = 8.0,
                    subfault_width: float = 8.0) -> dtopotools.Fault:
    """
    Subdivide each subfault into an nstrike x ndip grid, where:
      nstrike = int(length / subfault_length)
      ndip    = int(width  / subfault_width)
    Uses rake = strike - 60° - 180° (plate motion assumption).
    Returns a new Fault with the subdivided subfaults.
    """
    phi_plate = 60.0  # plate motion angle clockwise from north
    new_subfaults = []

    for sf in fault.subfaults:
        sf.rake = sf.strike - phi_plate - 180.0
        nstrike = int(sf.length / subfault_length)
        ndip = int(sf.width / subfault_width)

        if nstrike > 0 and ndip > 0:
            subdivided = dtopotools.SubdividedPlaneFault(sf, nstrike, ndip)
            new_subfaults.extend(subdivided.subfaults)
        else:
            # If the block is smaller than the target size, keep it as-is
            new_subfaults.append(sf)

    return dtopotools.Fault(subfaults=new_subfaults)

def generate_kl_scenarios(fault_geometry, n_scenarios=100, n_terms=20):
    """
    Generate scenarios using K-L expansion instead of k-means
    """
    # 1. Compute correlation matrix from fault geometry
    D, Dstrike, Ddip = compute_subfault_distances(fault_geometry)
    C = compute_correlation_matrix(Dstrike, Ddip, Lstrike, Ldip)
    
    # 2. Compute eigenmodes
    lam, V, mean_slip_g = compute_eigenmodes(C, mean_slip)
    
    # 3. Generate scenarios using K-L expansion
    scenarios = []
    for i in range(n_scenarios):
        z = np.random.randn(n_terms)
        slip = KL_expansion(z, lam, V, fault_geometry, tau)
        scenarios.append(slip)
    
    return scenarios


def main():
    print("=== KL-Expansion Test ===")

    # 1) Load CSZ fault geometry
    try:
        fault = setup_fault_from_csv('CSZe01.csv')
        print(f"Loaded fault with {len(fault.subfaults)} subfaults")
    except Exception as e:
        print(f"Error loading CSZe01.csv: {e}")
        return

    # 2) Subdivide fault (8 km x 8 km blocks, adjust if needed)
    subdivided_fault = subdivide_fault(fault, subfault_length=8.0, subfault_width=8.0)
    print(f"Subdivided fault has {len(subdivided_fault.subfaults)} subfaults")

    # 3) Generate a small batch of scenarios for a quick test
    #    (use small n_scenarios/n_terms to keep memory and runtime low)
    scenarios, weights = generate_kl_scenarios(
        subdivided_fault,
        n_scenarios=5,
        n_terms=10,
        Lstrike=130e3,
        Ldip=40e3,
        Mw_desired=9.0,
    )

    print(f"Generated {scenarios.shape[0]} scenarios with {scenarios.shape[1]} subfaults each")
    print(f"Slip range: {scenarios.min():.3f} to {scenarios.max():.3f} m")
    print(f"Weight sum: {weights.sum():.6f}")

    # 4) Quick plot: first scenario slip along subfault index
    plt.figure(figsize=(10, 4))
    plt.plot(scenarios[0], marker='o', ms=3)
    plt.title('KL Scenario #1 (slip per subfault)')
    plt.xlabel('Subfault Index')
    plt.ylabel('Slip (m)')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('kl_expansion_test_scenario1.png', dpi=150, bbox_inches='tight')
    print("Saved plot: kl_expansion_test_scenario1.png")

    # 5) Mean slip across scenarios (per subfault)
    mean_slip = scenarios.mean(axis=0)
    plt.figure(figsize=(10, 4))
    plt.plot(mean_slip, color='tab:red')
    plt.title('Mean Slip Across KL Scenarios')
    plt.xlabel('Subfault Index')
    plt.ylabel('Mean Slip (m)')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('kl_expansion_test_mean.png', dpi=150, bbox_inches='tight')
    print("Saved plot: kl_expansion_test_mean.png")

    print("=== KL-Expansion Test: DONE ===")

if __name__ == "__main__":
    main()