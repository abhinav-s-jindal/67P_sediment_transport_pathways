################################################################################
#                            Shape Model Files                                 #
################################################################################

Overview
--------
3D shape model of comet 67P/Churyumov-Gerasimenko used for sediment transport
simulations. These files define the surface geometry onto which ejection 
probabilities and particle distributions are mapped.

################################################################################
#                                   Files                                      #
################################################################################

CG_s7v1_10Km_16m_13M_B_INLET.obj
---------------------------------
Format: Wavefront OBJ (ASCII)

Description: High-resolution 3D mesh of comet 67P nucleus. Contains vertex 
positions (v), vertex normals (vn), and triangular faces (f) defining the 
surface geometry.

Model Statistics:
	Vertices: 220,300
	Faces: 440,596
	Resolution: ~10 meters per facet

Units: kilometers

Face Indexing: Faces are indexed sequentially starting from 1. 

Usage: 
	- Base geometry for all transport simulations
	- Visualization substrate for probability and distribution maps
	- Reference for facet-based analysis and source-target relationships

Compatibility: Can be opened in ParaView, MeshLab, Blender, or loaded 
programmatically using Python libraries (trimesh, pyvista, etc.)