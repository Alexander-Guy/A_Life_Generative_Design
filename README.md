Generative design with artifical life                                                          

Author        : Alexander Guy

Date          : 12/09/2020

Version       : 1.1.2

Software used : MATLAB 2019a, Solidworks 2018, 3D painter, 3D builder

Instructions :

Hello, welcome to the 3D artifical life based approach to mechanical component generation and optimisation, the directory 
is fairly straightfoward. Each Matlab code file exectutes either a biologically inspired algorithm or the STL exporter
function.The slime and swarm methods export a matrix file which the STL exported can then open and export as an STL file.

The leaf method has it's STL exporter tied into the main file, this is because the mesh size for the leaf is so small that
importing and exporting takes much longer. To operate the leaf generation code, please open in MATLAB 2019a or later and execute,
please ensure the geometry directory is in the location it was when downloaded. After the leaf algorithm has been completed an STL file
will be sitting in the same folder as the code file, please open this with any 3D printing or visualisation software (3D painter is a 
windows default). 

The slime mold or swarm code files need to be opened and run, after they are complete there will be a '.m' file ready for the STL exporter 
sitting in the smae directory as the code files. Please open the STL exporter and run it, follow instructions to export the resulting STL
file and open it with the afformentioned 3D visualise tools. 

If you wish to look at and examine the stress/3d recreated generated components please navigate to the 'geometry/recreated_parts' folder and 
open the '.SLDPT' files with Solidworks 2018 or later. On the feature tree to the left of the screen you will see a test result folder, click on this,
then the result types to see the heat map projections of stress/displacent/strain analysis. 
