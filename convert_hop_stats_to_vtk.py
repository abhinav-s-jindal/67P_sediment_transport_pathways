"""
VTK Converter for Particle Distribution Visualization
======================================================

This script converts particle distribution text files to VTK format for 3D
visualization in ParaView or other VTK-compatible software. Each facet's
particle count is mapped as a scalar field onto the comet surface geometry.

Required Input Files:
--------------------
1. Shape model (OBJ): CG_s7v1_10Km_16m_13M_B_INLET.obj
   - 3D mesh geometry with vertices and faces
   - Located in: obj_files/

2. Particle distribution files (TXT): particle_distribution_v0x_hopN.txt
   - Output from multi-hop redistribution simulation
   - Located in: particle_distribution_outputs/v0x/txts/

Output:
-------
VTK PolyData files (ASCII format) with particle counts as cell data:
  - particle_distribution_v0x_init.vtk (initial state)
  - particle_distribution_v0x_hop1.vtk through hopN.vtk
Saved to: particle_distribution_outputs/v0x/vtk/

Visualization:
-------------
Open VTK files in ParaView and color by 'particle_count' scalar field.
Facets with zero particles are excluded from the color scale for better
contrast of active redistribution regions.

Author: Abhinav S. Jindal
"""

import numpy as np
import os

# ============================================================================
# Configuration
# ============================================================================

SHAPE_MODEL_FILE = "obj_files/CG_s7v1_10Km_16m_13M_B_INLET.obj"
PARTICLE_DIST_BASE_DIR = "particle_distribution_outputs"

# ============================================================================
# Core Functions
# ============================================================================

def read_obj_file(obj_filename):
    """
    Read OBJ file and extract vertices and faces
    
    Parameters:
    -----------
    obj_filename : str
        Path to OBJ file
    
    Returns:
    --------
    vertices : numpy array
        Nx3 array of vertex coordinates
    faces : list
        List of face vertex indices (0-indexed)
    """
    vertices = []
    faces = []
    
    with open(obj_filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('v '):
                parts = line.split()
                vertex = [float(parts[1]), float(parts[2]), float(parts[3])]
                vertices.append(vertex)
            elif line.startswith('f '):
                parts = line.split()
                face_vertices = []
                for part in parts[1:]:
                    vertex_idx = int(part.split('/')[0]) - 1  # Convert to 0-indexed
                    face_vertices.append(vertex_idx)
                faces.append(face_vertices)
    
    return np.array(vertices), faces

def read_particle_distribution(txt_filename):
    """
    Read particle distribution from text file
    
    Parameters:
    -----------
    txt_filename : str
        Path to particle distribution text file
    
    Returns:
    --------
    numpy array : Particle counts per facet
    """
    particle_counts = []
    with open(txt_filename, 'r') as f:
        for line in f:
            count = float(line.strip())
            particle_counts.append(count)
    
    return np.array(particle_counts)

def write_vtk_file(vertices, faces, particle_data, vtk_filename, exclude_zeros=True):
    """
    Write VTK PolyData file with particle counts as cell data
    
    Parameters:
    -----------
    vertices : numpy array
        Vertex coordinates
    faces : list
        Face vertex indices
    particle_data : numpy array
        Particle count per facet
    vtk_filename : str
        Output VTK file path
    exclude_zeros : bool
        If True, zero-particle facets are excluded from color scale
    """
    if exclude_zeros:
        non_zero_mask = particle_data > 0
        filtered_data = particle_data[non_zero_mask]
        if len(filtered_data) == 0:
            print(f"  Warning: All particle values are zero")
            filtered_data = particle_data
    else:
        filtered_data = particle_data
    
    if exclude_zeros:
        num_excluded = np.sum(particle_data == 0)
        print(f"  Excluded {num_excluded} faces with 0.00 particles")
        print(f"  Color range: {np.min(filtered_data):.2f} to {np.max(filtered_data):.2f}")
    
    with open(vtk_filename, 'w') as f:
        # VTK header
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Particle Distribution Mesh\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        
        # Write vertices
        f.write(f"POINTS {len(vertices)} float\n")
        for vertex in vertices:
            f.write(f"{vertex[0]:.6f} {vertex[1]:.6f} {vertex[2]:.6f}\n")
        
        # Write faces
        total_face_entries = sum(len(face) + 1 for face in faces)
        f.write(f"POLYGONS {len(faces)} {total_face_entries}\n")
        for face in faces:
            f.write(f"{len(face)} " + " ".join(map(str, face)) + "\n")
        
        # Write particle data as cell data
        f.write(f"CELL_DATA {len(faces)}\n")
        f.write("SCALARS particle_count float 1\n")
        f.write("LOOKUP_TABLE default\n")
        
        for count in particle_data:
            if exclude_zeros and count == 0.0:
                f.write("0.0\n")
            else:
                f.write(f"{count:.6f}\n")

def convert_particle_files_to_vtk(velocity, obj_filename=SHAPE_MODEL_FILE, 
                                  base_dir=PARTICLE_DIST_BASE_DIR):
    """
    Convert all particle distribution files for a given velocity to VTK format
    
    Parameters:
    -----------
    velocity : float
        Ejection velocity in m/s (0.1 to 0.9)
    obj_filename : str
        Path to shape model OBJ file
    base_dir : str
        Base directory containing particle distribution outputs
    """
    velocity_str = f"{int(velocity * 10):02d}"
    velocity_label = f"v{velocity_str}"
    input_dir = os.path.join(base_dir, velocity_label)
    txt_dir = os.path.join(input_dir, "txts")
    vtk_output_dir = os.path.join(input_dir, "vtk")
    
    os.makedirs(vtk_output_dir, exist_ok=True)
    
    print(f"Converting particle distributions to VTK for velocity {velocity} m/s")
    print(f"Shape model: {obj_filename}")
    print(f"Input directory: {txt_dir}")
    print(f"Output directory: {vtk_output_dir}\n")
    
    # Read shape model
    try:
        vertices, faces = read_obj_file(obj_filename)
        print(f"Loaded mesh: {len(vertices)} vertices, {len(faces)} faces\n")
    except Exception as e:
        print(f"Error reading OBJ file: {e}")
        return
    
    # Find all particle distribution files
    files_to_process = []
    
    init_file = os.path.join(txt_dir, f"particle_distribution_{velocity_label}_init.txt")
    if os.path.exists(init_file):
        files_to_process.append(("init", init_file))
    
    for hop in range(1, 1001):
        hop_file = os.path.join(txt_dir, f"particle_distribution_{velocity_label}_hop{hop}.txt")
        if os.path.exists(hop_file):
            files_to_process.append((f"hop{hop}", hop_file))
    
    if not files_to_process:
        print("No particle distribution files found!")
        return
    
    print(f"Found {len(files_to_process)} files to convert\n")
    
    # Convert each file
    for label, txt_file in files_to_process:
        print(f"Processing {label}...")
        
        try:
            particle_data = read_particle_distribution(txt_file)
            
            # Verify data length matches face count
            if len(particle_data) != len(faces):
                print(f"  Warning: Data length mismatch ({len(particle_data)} vs {len(faces)} faces)")
                if len(particle_data) > len(faces):
                    particle_data = particle_data[:len(faces)]
                else:
                    padded_data = np.zeros(len(faces))
                    padded_data[:len(particle_data)] = particle_data
                    particle_data = padded_data
            
            vtk_filename = os.path.join(vtk_output_dir, 
                                       f"particle_distribution_{velocity_label}_{label}.vtk")
            
            write_vtk_file(vertices, faces, particle_data, vtk_filename, exclude_zeros=True)
            print(f"  Saved: {vtk_filename}\n")
            
        except Exception as e:
            print(f"  Error processing {label}: {e}\n")
    
    print(f"VTK conversion complete!")
    print(f"Files saved in: {vtk_output_dir}")
    print("Load in ParaView and color by 'particle_count' scalar field")

# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    velocity = float(input("Enter ejection velocity (0.1-0.9 m/s): "))
    convert_particle_files_to_vtk(velocity)
