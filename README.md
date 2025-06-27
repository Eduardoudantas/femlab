# Matlab-femlab

## Overview

This project is a didactic program for **Linear and Dynamic Analysis using the Finite Element Method (FEM)**, developed for educational purposes in Civil Engineering. It is designed to analyze the static and dynamic behavior of plate structures using quadrilateral finite elements (Q4, Q8, Q9) with Mindlin-Reissner theory.

- **Author:** Eduardo Uchôa Dantas
- **Institution:** Universidade Federal do Pará (UFPa), ITEC

## Features

- Static and modal (vibration) analysis of plate structures.
- Support for Q4, Q8, and Q9 quadrilateral elements.
- Customizable mesh, material properties, loads, and boundary conditions.
- Visualization of undeformed/deformed structures and modal shapes.

## Project Structure

- `AnalisePrincipalFEM.m`: **Main analysis function**. Assembles the global system, solves for displacements, reactions, stresses, strains, and natural frequencies.
- `Inputq9_Model.m`, `Inputq9_Mierovich.m`, `Inputq9_Menao.m`, `Inputq4model.m`: **Input scripts**. Define geometry, mesh, material properties, loads, and call the main analysis.
- `ElemQuadrilateralq4mindlin.m`, `ElemQuadrilateralq8mindlin.m`, `ElemQuadrilateralq9mindlin.m`: Element routines for Q4, Q8, and Q9 elements (stiffness, mass, loads, stresses, strains).
- `PtsGauss1d.m`, `PtsGauss2d.m`: Gauss integration points for numerical integration.
- `Matrizjacobiana2d.m`: Jacobian matrix calculation for coordinate transformation.
- `desenhaestinderformadaq93d.m`, `desenhaesdeformadaq93drev.m`, `desenhaformasmodais.m`: Visualization scripts for mesh, deformed shape, and modal shapes.

## How the Program Works

1. **Define the Model**: Use one of the input scripts (e.g., `Inputq9_Model.m`) to set up the plate dimensions, mesh density, material properties, boundary conditions, and loads.
2. **Mesh Generation**: The input script generates node coordinates and element connectivity automatically.
3. **Boundary Conditions & Loads**: Specify which nodes are fixed or loaded.
4. **Element Data**: Assign material and geometric properties to each element.
5. **Run Analysis**: The input script calls `AnalisePrincipalFEM`, which:
   - Assembles global stiffness and mass matrices.
   - Applies boundary conditions and loads.
   - Solves for nodal displacements and reactions.
   - Computes element stresses and strains.
   - Performs modal analysis (natural frequencies and mode shapes).
6. **Results & Visualization**: Results are displayed in the MATLAB console and visualized using the provided plotting scripts.

## How to Use

1. **Open MATLAB** and navigate to the project folder.
2. **Choose and edit an input script** (e.g., `Inputq9_Model.m`) to match your problem setup.
3. **Run the input script** in the MATLAB command window:
   ```
   Inputq9_Model
   ```
   or for Q4 elements:
   ```
   Inputq4model
   ```
4. **View results** in the MATLAB console and graphical figures.

## Customization

- **Mesh Density**: Change `NumElementosx` and `NumElementosy` in the input script.
- **Material Properties**: Edit `E`, `nu`, `t`, etc.
- **Loads**: Modify the `cargasNos` matrix for different load cases.
- **Boundary Conditions**: Adjust the `restrs` matrix to fix or release nodes as needed.
- **Element Type**: Set `tipoElem` to `'Quadrilateralq4mindlin'`, `'Quadrilateralq8mindlin'`, or `'Quadrilateralq9mindlin'`.

## Example

To analyze a Q9 plate, open `Inputq9_Model.m`, adjust the parameters, and run:
```matlab
Inputq9_Model
```
You will see:
- Node coordinates and connectivity
- Nodal displacements and reactions
- Element stresses and strains
- Natural frequencies
- Plots of the undeformed/deformed structure and modal shapes

## Requirements

- MATLAB (no toolboxes required)
- Basic knowledge of FEM and MATLAB scripting

## Educational Notes

- The code is heavily commented for learning purposes.
- Each function and script is modular and can be studied independently.
- The project is suitable for students and educators in structural analysis and FEM courses. 