"""
Cleanup procedures to make last mods to par files
"""
import os

fid = "EXPORT/meshfem3D_files/Mesh_Par_file"
with open(fid, "r") as f:
    lines = f.readlines()

with open(fid, "w") as f:
    for line in lines[:-1]:
        if line.startswith("NMATERIALS"):
            new_line = line.replace("1", "2")
            f.write(new_line)
            f.write("-1 tomography elastic tomography_model_shallow.xyz 0 2\n")
        elif line.startswith("-1 tomography"):
            f.write("-2 tomography elastic tomography_model_whole.xyz 0 2\n")
        elif line.startswith("1       "):
            new_line = line.replace("-1", "-2")
            f.write(line.replace("-1", "-2"))
        elif line.startswith("# ADD"):
            f.write(line[2:])
        elif line.startswith("# NUMBER_"):
            f.write(line[2:])
        else:
            f.write(line)
    f.write(lines[-1])




            
