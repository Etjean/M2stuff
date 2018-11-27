#!/bin/bash

# Remove waters from the PDB
grep -v HOH 1aki.pdb > 1AKI_clean.pdb

# Choose the force field. 
# GRO file : structure in the force field (atom positions)
# TOP file : topology for the force field (atom & bonds characteristics)
# ITP file : positions restraints of the non-hydrogen atoms
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce

# Define box dimensions
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic

# Fill the box with water
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

# Add counter-ions
# Generate ions.tpr, binary file from the parameters in ions.mdp
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Minimization
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
# analysis
gmx energy -f em.edr -o potential.xvg






