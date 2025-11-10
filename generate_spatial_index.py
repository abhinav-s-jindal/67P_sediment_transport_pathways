"""
Spatial Index Generator for Comet 67P Shape Model
==================================================

This script generates a KDTree spatial index for efficient facet lookup operations
in the interactive reverse trajectory visualizer. The spatial index enables fast
neighbor queries for the brush selection tool.

Required Input Files:
--------------------
Shape model (OBJ): CG_s7v1_10Km_16m_13M_B_INLET.obj
   - 3D mesh geometry
   - Located in: obj_files/

Output:
-------
cell_spatial_index.pkl: Pickle file containing cell centers and KDTree
   - Used by trajectory_source_visualizer.py
   - Enables fast facet neighbor queries

Usage:
------
Run this script once before using the interactive visualizer:
    python generate_spatial_index.py

Requirements:
-------------
- VTK (for OBJ loading)
- NumPy (for array operations)
- scipy (for KDTree construction)

Author: Abhinav S. Jindal
"""

import vtk
import numpy as np
import pickle
from scipy.spatial import cKDTree
import time
import os

# ============================================================================
# Configuration
# ============================================================================

SHAPE_MODEL_FILE = "obj_files/CG_s7v1_10Km_16m_13M_B_INLET.obj"
OUTPUT_FILE = "cell_spatial_index.pkl"

# ============================================================================
# Core Functions
# ============================================================================

def load_shape_model(obj_filename):
    """
    Load comet shape model from OBJ file
    
    Parameters:
    -----------
    obj_filename : str
        Path to OBJ file
        
    Returns:
    --------
    vtkPolyData : Loaded mesh data
    """
    print(f"Loading shape model from {obj_filename}...")
    
    reader = vtk.vtkOBJReader()
    reader.SetFileName(obj_filename)
    reader.Update()
    
    polydata = reader.GetOutput()
    num_cells = polydata.GetNumberOfCells()
    num_points = polydata.GetNumberOfPoints()
    
    print(f"Loaded: {num_points:,} vertices, {num_cells:,} facets\n")
    
    return polydata

def extract_cell_centers(polydata):
    """
    Extract centroid coordinates for all mesh facets
    
    Parameters:
    -----------
    polydata : vtkPolyData
        Input mesh data
        
    Returns:
    --------
    numpy.ndarray : Nx3 array of cell center coordinates
    """
    print("Extracting cell centers...")
    
    num_cells = polydata.GetNumberOfCells()
    centers = np.zeros((num_cells, 3))
    
    report_interval = max(1, num_cells // 20)
    
    for i in range(num_cells):
        if i % report_interval == 0:
            progress = (i / num_cells) * 100
            print(f"  Progress: {progress:.1f}%")
        
        cell = polydata.GetCell(i)
        points = cell.GetPoints()
        num_points = points.GetNumberOfPoints()
        
        center = [0.0, 0.0, 0.0]
        for j in range(num_points):
            point = points.GetPoint(j)
            for k in range(3):
                center[k] += point[k]
        
        for k in range(3):
            center[k] /= num_points
            
        centers[i] = center
    
    print(f"Extracted {num_cells:,} cell centers\n")
    return centers

def create_spatial_tree(centers):
    """
    Build KDTree spatial index for fast neighbor queries
    
    Parameters:
    -----------
    centers : numpy.ndarray
        Array of cell center coordinates
        
    Returns:
    --------
    scipy.spatial.cKDTree : Spatial tree for neighbor queries
    """
    print("Building KDTree spatial index...")
    start_time = time.time()
    
    tree = cKDTree(centers)
    
    elapsed = time.time() - start_time
    print(f"KDTree built in {elapsed:.2f}s")
    print(f"Indexed {len(centers):,} points\n")
    
    return tree

def save_spatial_index(centers, tree, output_filename):
    """
    Save spatial index to pickle file
    
    Parameters:
    -----------
    centers : numpy.ndarray
        Cell center coordinates
    tree : cKDTree
        Spatial tree
    output_filename : str
        Output pickle file path
    """
    print(f"Saving spatial index to {output_filename}...")
    
    spatial_data = {
        'centers': centers,
        'tree': tree,
        'num_facets': len(centers),
        'creation_time': time.strftime('%Y-%m-%d %H:%M:%S')
    }
    
    with open(output_filename, 'wb') as f:
        pickle.dump(spatial_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    file_size = get_file_size(output_filename)
    print(f"Saved successfully ({file_size})\n")

def get_file_size(filename):
    """Convert file size to human-readable format."""
    size_bytes = os.path.getsize(filename)
    
    if size_bytes < 1024:
        return f"{size_bytes} bytes"
    elif size_bytes < 1024**2:
        return f"{size_bytes/1024:.1f} KB"
    elif size_bytes < 1024**3:
        return f"{size_bytes/(1024**2):.1f} MB"
    else:
        return f"{size_bytes/(1024**3):.1f} GB"

def test_spatial_index(centers, tree):
    """
    Test spatial index with sample queries
    
    Parameters:
    -----------
    centers : numpy.ndarray
        Cell center coordinates  
    tree : cKDTree
        Spatial tree
    """
    print("Testing spatial index...")
    
    test_center = centers[0]
    test_radius = 0.1
    
    start_time = time.time()
    nearby_indices = tree.query_ball_point(test_center, test_radius)
    query_time = time.time() - start_time
    
    print(f"  Test query at facet 0:")
    print(f"    Radius: {test_radius}")
    print(f"    Found: {len(nearby_indices)} neighbors")
    print(f"    Query time: {query_time*1000:.2f} ms")
    
    import random
    for i in range(3):
        rand_idx = random.randint(0, len(centers)-1)
        test_center = centers[rand_idx]
        nearby_indices = tree.query_ball_point(test_center, test_radius)
        print(f"  Test {i+1}: facet {rand_idx} has {len(nearby_indices)} neighbors")
    
    print()

# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Generate spatial index for comet shape model."""
    
    print("=" * 70)
    print("COMET 67P SPATIAL INDEX GENERATOR")
    print("=" * 70)
    print(f"Input:  {SHAPE_MODEL_FILE}")
    print(f"Output: {OUTPUT_FILE}\n")
    
    start_time = time.time()
    
    try:
        polydata = load_shape_model(SHAPE_MODEL_FILE)
        centers = extract_cell_centers(polydata)
        tree = create_spatial_tree(centers)
        save_spatial_index(centers, tree, OUTPUT_FILE)
        test_spatial_index(centers, tree)
        
        total_time = time.time() - start_time
        print("=" * 70)
        print("SPATIAL INDEX GENERATION COMPLETE")
        print("=" * 70)
        print(f"Facets processed: {len(centers):,}")
        print(f"Total time: {total_time:.2f}s")
        print(f"Output file: {OUTPUT_FILE}\n")
        print("Ready to run trajectory_source_visualizer.py!")
        
    except FileNotFoundError:
        print(f"\nERROR: Could not find {SHAPE_MODEL_FILE}")
        print("Ensure the shape model is in the obj_files/ directory")
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()