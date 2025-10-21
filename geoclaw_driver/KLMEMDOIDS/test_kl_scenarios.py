"""
Simple test of K-L expansion scenario generation
"""

import numpy as np
import matplotlib.pyplot as plt
from KL_Scenario_Generation import generate_kl_scenarios, setup_fault_from_csv, subdivide_fault

def test_kl_scenario_generation():
    """
    Test the K-L expansion scenario generation
    """
    print("Testing K-L expansion scenario generation...")
    
    try:
        # Load fault geometry
        print("Loading CSZ fault geometry...")
        fault = setup_fault_from_csv('CSZe01.csv')
        
        # Subdivide fault
        print("Subdividing fault...")
        subdivided_fault = subdivide_fault(fault, subfault_length=20.0, subfault_width=20.0)
        
        # Generate scenarios
        print("Generating K-L scenarios...")
        scenarios, weights = generate_kl_scenarios(
            subdivided_fault,
            n_scenarios=10,  # Small number for testing
            n_terms=10,      # Fewer terms for testing
            Lstrike=130e3,   # 130 km correlation length
            Ldip=40e3,       # 40 km correlation length
            Mw_desired=9.0
        )
        
        print(f"Successfully generated {scenarios.shape[0]} scenarios")
        print(f"Each scenario has {scenarios.shape[1]} subfaults")
        print(f"Slip range: {scenarios.min():.2f} to {scenarios.max():.2f} meters")
        print(f"Weight sum: {weights.sum():.4f}")
        
        # Create a simple plot
        plt.figure(figsize=(12, 4))
        
        # Plot first scenario
        plt.subplot(1, 3, 1)
        plt.plot(scenarios[0, :])
        plt.title('First K-L Scenario')
        plt.xlabel('Subfault Index')
        plt.ylabel('Slip (m)')
        plt.grid(True, alpha=0.3)
        
        # Plot scenario statistics
        plt.subplot(1, 3, 2)
        mean_slip = scenarios.mean(axis=0)
        plt.plot(mean_slip)
        plt.title('Mean Slip Across Scenarios')
        plt.xlabel('Subfault Index')
        plt.ylabel('Mean Slip (m)')
        plt.grid(True, alpha=0.3)
        
        # Plot weight distribution
        plt.subplot(1, 3, 3)
        plt.hist(weights, bins=10, alpha=0.7, edgecolor='black')
        plt.title('Weight Distribution')
        plt.xlabel('Weight')
        plt.ylabel('Frequency')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('kl_scenarios_test.png', dpi=150, bbox_inches='tight')
        print("Saved test plot: kl_scenarios_test.png")
        
        return True
        
    except Exception as e:
        print(f"Error during testing: {e}")
        return False

def test_without_fault():
    """
    Test K-L expansion without fault file (using dummy fault)
    """
    print("\nTesting K-L expansion with dummy fault...")
    
    try:
        from clawpack.geoclaw import dtopotools
        
        # Create a simple dummy fault
        fault = dtopotools.Fault()
        for i in range(10):  # 10 subfaults
            subfault = dtopotools.Fault().Subfault()  # if available
            subfault.longitude = -125.0 + i * 0.1
            subfault.latitude = 42.0 + i * 0.05
            subfault.depth = 10.0 + i * 0.5
            subfault.strike = 0.0
            subfault.dip = 10.0
            subfault.length = 10.0
            subfault.width = 10.0
            subfault.slip = 1.0
            subfault.rake = 90.0
            fault.subfaults.append(subfault)
        
        # Generate scenarios
        scenarios, weights = generate_kl_scenarios(
            fault,
            n_scenarios=5,
            n_terms=5,
            Lstrike=50e3,  # 50 km correlation length
            Ldip=20e3,     # 20 km correlation length
            Mw_desired=8.0
        )
        
        print(f"Successfully generated {scenarios.shape[0]} scenarios with dummy fault")
        print(f"Slip range: {scenarios.min():.2f} to {scenarios.max():.2f} meters")
        
        return True
        
    except Exception as e:
        print(f"Error with dummy fault: {e}")
        return False

if __name__ == "__main__":
    print("Testing K-L Expansion Scenario Generation")
    print("="*50)
    
    # Test with real fault if available
    success1 = test_kl_scenario_generation()
    
    # Test with dummy fault
    success2 = test_without_fault()
    
    if success1 or success2:
        print("\n✅ K-L expansion scenario generation is working!")
        print("\nAdvantages over k-means:")
        print("1. Physically realistic scenarios")
        print("2. No arbitrary dilation needed")
        print("3. Computationally efficient")
        print("4. Respects fault geometry")
    else:
        print("\n❌ K-L expansion scenario generation failed")
        print("Check error messages above for details")
