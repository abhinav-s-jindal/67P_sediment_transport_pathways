"""
Multi-Hop Sediment Redistribution Simulation for Comet 67P
===========================================================

This script simulates iterative sediment redistribution across the comet surface
through ballistic particle transport. Starting from an initial distribution,
particles are ejected from facets according to ejection probabilities and land
on destination facets following pre-computed transport pathways.

Required Input Files:
--------------------
1. Transport dictionary (PKL): transport_dict_velocity_0x_shapemodel.pkl
   - Contains source-to-destination facet mappings for a given velocity
   - Located in: transport_dict_files/
   
2. Ejection probability map (TXT): CG_s7v1_10Km_16m_13M_B_INLET_ejection_probability.txt
   - Normalized ejection probability for each facet (one value per line)
   - Located in: ejection_probability_maps/
   
3. Shape model (OBJ): CG_s7v1_10Km_16m_13M_B_INLET.obj
   - Reference geometry with 440,596 faces
   - Located in: obj_files/

Algorithm:
----------
Uses a two-phase processing approach to ensure particle conservation:
  Phase 1 (Read): Calculate all ejections from current distribution
  Phase 2 (Write): Apply all transfers simultaneously

This prevents order-dependency issues where facets that eject early in an
iteration would have their newly-received particles incorrectly re-ejected
in the same iteration.

Output:
-------
Creates text files for each hop iteration:
  - particle_distribution_v0x_init.txt (initial state)
  - particle_distribution_v0x_hop1.txt through hopN.txt
Each file contains 440,596 lines where line number = facet ID particle count.

Author: Abhinav S. Jindal
"""

import pickle
import numpy as np
from collections import defaultdict, Counter
import os

# ============================================================================
# Configuration
# ============================================================================

# Directory structure
TRANSPORT_DICT_DIR = "transport_dict_files"
EJECTION_PROB_DIR = "ejection_probability_maps"
SHAPE_MODEL_FILE = "obj_files/CG_s7v1_10Km_16m_13M_B_INLET.obj"

# Shape model parameters
MAX_FACET_ID = 440595  # 440,596 faces (0-indexed)

# ============================================================================
# Core Functions
# ============================================================================

def load_transport_data(velocity):
    """
    Load transport dictionary for specified velocity
    
    Parameters:
    -----------
    velocity : float
        Ejection velocity in m/s (0.1 to 0.9)
    
    Returns:
    --------
    dict : Transport dictionary mapping source facets to destinations
    """
    velocity_str = f"{int(velocity * 10):02d}"
    pickle_file = os.path.join(TRANSPORT_DICT_DIR, 
                               f"transport_dict_velocity_{velocity_str}_shapemodel.pkl")
    
    with open(pickle_file, 'rb') as f:
        transport_dict = pickle.load(f)
    return transport_dict

def load_ejection_probabilities(prob_file):
    """
    Load ejection probabilities from text file
    
    Parameters:
    -----------
    prob_file : str
        Path to ejection probability file
    
    Returns:
    --------
    dict : Facet ID (0-indexed) to ejection probability mapping
    """
    ejection_probs = {}
    with open(prob_file, 'r') as f:
        for line_num, line in enumerate(f):
            prob = float(line.strip())
            ejection_probs[line_num] = prob
    return ejection_probs

def initialize_particle_distribution(transport_dict):
    """
    Initialize particle distribution based on connectivity.
    Each facet receives particles equal to its total outgoing connections.
    Missing facets are padded with zeros.
    
    Parameters:
    -----------
    transport_dict : dict
        Transport dictionary from load_transport_data()
    
    Returns:
    --------
    dict : Facet ID to particle count mapping
    """
    particle_distribution = {facet_id: 0.0 for facet_id in range(MAX_FACET_ID + 1)}
    
    for source_facet, destinations_dict in transport_dict.items():
        if isinstance(destinations_dict, dict) and 'destination_facets' in destinations_dict:
            destinations = destinations_dict['destination_facets']
        else:
            destinations = destinations_dict
        
        particle_distribution[source_facet] = float(len(destinations))
    
    return particle_distribution

def perform_single_hop_two_phase(current_distribution, transport_dict, ejection_probs):
    """
    Perform one redistribution iteration using two-phase approach
    
    Phase 1: Calculate all ejections (read-only)
    Phase 2: Apply all transfers simultaneously (write)
    
    This ensures particle conservation and eliminates order-dependency.
    
    Parameters:
    -----------
    current_distribution : dict
        Current particle counts per facet
    transport_dict : dict
        Transport pathways
    ejection_probs : dict
        Ejection probabilities per facet
    
    Returns:
    --------
    dict : Updated particle distribution after one hop
    """
    ejection_matrix = defaultdict(lambda: defaultdict(float))
    
    # Phase 1: Calculate all ejections
    for source_facet, current_particles in current_distribution.items():
        if current_particles <= 0 or source_facet not in transport_dict:
            continue
        
        prob_ejection = ejection_probs.get(source_facet, 0.0)
        if prob_ejection <= 0:
            continue
        
        destinations_dict = transport_dict[source_facet]
        if isinstance(destinations_dict, dict) and 'destination_facets' in destinations_dict:
            destinations = destinations_dict['destination_facets']
        else:
            destinations = destinations_dict
        
        total_destination_count = len(destinations)
        destination_counts = Counter(destinations)
        particles_available_for_ejection = current_particles * prob_ejection
        
        for dest_facet, occurrences in destination_counts.items():
            particles_to_dest = (particles_available_for_ejection / total_destination_count) * occurrences
            ejection_matrix[source_facet][dest_facet] += particles_to_dest
    
    # Phase 2: Apply all changes
    new_distribution = current_distribution.copy()
    
    for source_facet, destinations in ejection_matrix.items():
        total_ejected = sum(destinations.values())
        new_distribution[source_facet] -= total_ejected
        
        for dest_facet, particles_received in destinations.items():
            if dest_facet not in new_distribution:
                new_distribution[dest_facet] = 0.0
            new_distribution[dest_facet] += particles_received
    
    return new_distribution

def save_distribution(distribution, filename):
    """
    Save particle distribution to text file (one value per line)
    Line number corresponds to facet ID
    
    Parameters:
    -----------
    distribution : dict
        Facet ID to particle count mapping
    filename : str
        Output file path
    """
    with open(filename, 'w') as f:
        for facet_id in range(MAX_FACET_ID + 1):
            if facet_id in distribution:
                f.write(f"{distribution[facet_id]:.2f}\n")
            else:
                f.write("0.00\n")

def run_particle_hop_analysis(velocity, ejection_prob_file, number_of_hops, output_base_dir):
    """
    Run complete multi-hop sediment redistribution simulation
    
    Parameters:
    -----------
    velocity : float
        Ejection velocity in m/s (0.1 to 0.9)
    ejection_prob_file : str
        Path to ejection probability file
    number_of_hops : int
        Number of redistribution iterations to perform
    output_base_dir : str
        Base directory for output files
    
    Returns:
    --------
    dict : Final particle distribution after all hops
    """
    if not (0.1 <= velocity <= 0.9):
        raise ValueError("Velocity must be between 0.1 and 0.9 m/s")
    
    velocity_str = f"{int(velocity * 10):02d}"
    velocity_label = f"v{velocity_str}"
    
    print(f"Running multi-hop analysis for velocity {velocity} m/s")
    print(f"Number of hops: {number_of_hops}")
    print("Using two-phase processing for particle conservation\n")
    
    output_dir = os.path.join(output_base_dir, velocity_label, 'txts')
    os.makedirs(output_dir, exist_ok=True)
    
    transport_dict = load_transport_data(velocity)
    ejection_probs = load_ejection_probabilities(ejection_prob_file)
    current_distribution = initialize_particle_distribution(transport_dict)
    
    # Save initial distribution
    init_file = os.path.join(output_dir, f"particle_distribution_{velocity_label}_init.txt")
    save_distribution(current_distribution, init_file)
    
    total_particles = sum(current_distribution.values())
    num_active_facets = sum(1 for count in current_distribution.values() if count > 0)
    print(f"Initial: {total_particles:.2f} particles across {num_active_facets} facets\n")
    
    # Perform hops
    for hop in range(1, number_of_hops + 1):
        print(f"Hop {hop}...", end=" ")
        current_distribution = perform_single_hop_two_phase(current_distribution, 
                                                            transport_dict, 
                                                            ejection_probs)
        
        hop_file = os.path.join(output_dir, f"particle_distribution_{velocity_label}_hop{hop}.txt")
        save_distribution(current_distribution, hop_file)
        
        total_after = sum(current_distribution.values())
        conservation = (total_after / total_particles) * 100
        print(f"Total: {total_after:.2f} | Conservation: {conservation:.4f}%")
    
    print(f"\nComplete! Output directory: {output_dir}")
    
    return current_distribution

# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    ejection_prob_file = os.path.join(EJECTION_PROB_DIR, 
                                     "CG_s7v1_10Km_16m_13M_B_INLET_ejection_probability.txt")
    
    output_base_dir = "particle_distribution_outputs"
    
    velocity = float(input("Enter ejection velocity (0.1-0.9 m/s): "))
    number_of_hops = int(input("Enter number of hops: "))
    
    final_distribution = run_particle_hop_analysis(
        velocity=velocity,
        ejection_prob_file=ejection_prob_file,
        number_of_hops=number_of_hops,
        output_base_dir=output_base_dir
    )
    
    print(f"\nOutput files: particle_distribution_v{int(velocity*10):02d}_init.txt through hop{number_of_hops}.txt")
    print("Each file has 440,596 lines where line N = facet N particle count")