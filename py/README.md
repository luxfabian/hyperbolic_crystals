# Python scripts

The scripts contained in this folder can be used to study the spectral properties of various hyperbolic lattices based on an existing `.words` file and the corresponding files for the regular representation of the group generators. We now give a brief overview of the files in this folder

## `iomodule.py`

Extracts the information from `group_specs.inp` and collects the data into a dictionary `access_point` which can be used by other scripts as a gateway to load the `.words` and/or the regular representation. It further contains the routine `fortan_format` which can be used to convert a `float` into a Fortran-readable string of the pre-defined E15.6 format.