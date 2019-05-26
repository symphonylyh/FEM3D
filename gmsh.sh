#!/bin/bash

# Gmsh is an open-sourced automatic 2D/3D mesh generation software: http://gmsh.info
# Some tutorials on Gmsh: http://gmsh.info/doc/texinfo/gmsh.html
# [1] How to generate structured grid (i.e. Q8 or B20 quad-like elements)? Keyword: Transfinite
#     2D: https://www.youtube.com/watch?v=O1FyiBBuN98
#     3D: https://www.youtube.com/watch?v=ewp3VGyymK4
# [2] How to refine mesh?
#     https://www.youtube.com/watch?v=kpWVNNHHdd8

# Step-by-step for the UI: "e" for end current operation, "q" for finish current function
# [1] Draw 
#     Geometry--Elementary entities--Add--Point-->Line-->Plane Surface--(Extrude)-->Volume
#     You can also add Physical groups to better organize things
# [2] Mesh 
#     Mesh--Define--Transfinite (Transfinite will make the mesh into quads)
#         Curve--Number of points (e.g. 3 points means 2 divisions)--Select lines
#         Surface--Select faces
#         Volume--Select centroid--Select vertices (NOTE: the node ordering MATTERs! you should follow a certain order for picking the nodes)
#     Mesh--Define--Recombine--Select surfaces (recombine will force the elements to be more structured)
#     Mesh--3D
# [3] Export
#     File--Export--Mesh-VTK--.vtk 
# As you can see, the above operations are just a few command lines in .geo file.
# We can also use script to further automate the generation by: http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options
export PATH=$PATH:/Applications/Gmsh.app/Contents/MacOS/
gmsh trial.geo -3 -format vtk -o trial.vtk # -3: 3D meshing
# .geo can also contain meshing commands, such as 'Mesh 3; Mesh.SecondOrderIncomplete=0;' Or pass as string '-string "Mesh.SecondOrderIncomplete=0;"' in command line
# More meshing commands: http://gmsh.info/doc/texinfo/gmsh.html#Mesh-options-list

# Be careful about the node ordering: http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering


# .geo file for the script 

# .vtk file for Paraview style output (.msh format is not human-readable)