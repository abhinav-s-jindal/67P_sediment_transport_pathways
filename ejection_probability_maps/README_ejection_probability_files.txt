################################################################################
#                        Ejection Probability Maps                             #
################################################################################

Overview
--------
Ejection probability distributions derived from the cumulative erosion map 
published in Groussin et al. (2025). These maps represent the normalized 
probability of particle ejection from each facet of the 67P/Churyumov-
Gerasimenko shape model.

Reference
---------
Erosion data source: Groussin et al. (2025)  
DOI: https://doi.org/10.1051/0004-6361/202452260

Shape Model
-----------
These probability maps correspond to the shape model:  
	../obj_files/CG_s7v1_10Km_16m_13M_B_INLET.obj

################################################################################
#                                   Files                                      #
################################################################################

CG_s7v1_10Km_16m_13M_B_INLET_ejection_probability.txt
------------------------------------------------------
Format: Plain text, one value per line  

Description: Normalized ejection probability for each facet of the shape model. 
Values represent the relative likelihood of particle ejection, derived from 
erosion rates. Facet indexing corresponds to the face ordering in the OBJ file.

Usage: Input for multi-hop sediment transport simulations.


CG_s7v1_10Km_16m_13M_B_INLET_ejection_probability.vtk
------------------------------------------------------
Format: VTK PolyData (ASCII)  

Description: 3D visualization file containing the same ejection probability 
data mapped onto the 67P surface geometry. Can be opened in ParaView or other 
VTK-compatible viewers.

Usage: Visualization and validation of ejection probability distribution across 
the comet surface.