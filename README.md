This project contains a collection of MATLAB scripts and associated data files used for performing various simulations and analyses related to farfield focusing, delta scaling tests, metaunit mapping, and phase plotting. Below is a brief description of each file and its purpose.

File Descriptions
Scripts
DeltaScalingTest.m

Description: This script performs delta scaling tests, which likely involve adjusting a scaling parameter (delta) and observing its effects on the simulation or model.
Usage: Run this script to perform the delta scaling test and generate relevant plots or output data.
farfieldFocusing.m

Description: This script is used to simulate and analyze farfield focusing, which is a technique used in optics and wave propagation.
Usage: Execute this script to perform farfield focusing simulations and visualize the results.
MetaunitMapper.m

Description: This script maps the properties of metaunits, which are the basic building blocks in metamaterials.
Usage: Use this script to map and analyze different metaunit configurations and their properties.
MetaunitMapperSquareStructure.m

Description: A variant of the MetaunitMapper script, tailored for square-structured metaunits.
Usage: Run this script to specifically analyze square metaunit structures.
nestedsweepD1D2_amplitude.m

Description: This script performs nested sweeps for two parameters (D1 and D2) and analyzes their effects on the amplitude of the system.
Usage: Execute this script to conduct the nested sweeps and plot the amplitude results.
nestedsweepD1D2_phaseplot.m

Description: Similar to the above, but focuses on generating phase plots instead of amplitude plots.
Usage: Use this script to generate phase plots from the nested sweeps of D1 and D2.
NetFarfield.m

Description: This script aggregates the results from various farfield simulations and provides a comprehensive analysis.
Usage: Run this script to compile and analyze the results from multiple farfield focusing experiments.
Data Files
deltacolormap.mat

Description: Contains colormap data used for visualizations in delta scaling tests.
Usage: Load this file in MATLAB to apply the colormap to your plots.
phasecolormap.mat

Description: Contains colormap data used for phase plotting.
Usage: Load this file in MATLAB to apply the colormap to your phase plots.
Autosave Files
farfieldFocusing.asv
Description: An autosave file associated with the farfieldFocusing.m script. Contains intermediate data or state information.
Usage: Typically used for recovery purposes. Load this file if you need to restore a previous session of farfieldFocusing.m.
Usage Instructions
Set Up MATLAB Environment:

Ensure MATLAB is installed and properly configured on your system.
Add the project directory to the MATLAB path using addpath('path_to_project_directory').
Running Scripts:

Open MATLAB and navigate to the project directory.
Run the desired script by typing its name in the command window (e.g., DeltaScalingTest).
Loading Data Files:

Use the load function to load .mat files as needed within your scripts.
Example: load('deltacolormap.mat').
Dependencies
MATLAB (R2021a or later recommended)
Signal Processing Toolbox (for some scripts)
Author
Name: [Your Name]
Contact: [Your Email]
License
This project is licensed under the MIT License. See the LICENSE file for details.

Acknowledgments
Acknowledge any individuals or organizations that contributed to the project.
