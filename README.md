# bachelor-thesis-discrete-minimal-surfaces

This repository contains a small collection of self-contained Julia scripts 
used to generate figures and exploratory visualizations 
with GLMakie for my Bachelor thesis on discrete minimal surfaces.

The code accompanies the thesis but is not required to understand the
mathematical results. All definitions, proofs, and statements in the thesis
are independent of the computational experiments.

## Contents

- `discrete_isothermic.jl`  
  Planar discrete conformal maps built by cr = -1 evolution starting from initial data prescribed on boundary vertices.
  Discrete I-minimal surfaces constructed via I-minimal Weierstrass representation applied to discrete conformal map.

- `discrete_s_isothermic.jl`  
  Planar Schramm circle pattern built by iteratively constructing tangent / orthogonal circle from initial data prescribed on the main diagonal.
  Discrete S-minimal surfaces visualized by spheres / circles constructed via S-minimal Weierstrass representation applied to Schramm circle pattern.

- `planar_variation.jl`  
  Interactive script visualizing behavior of planar discrete conformal map with half the vertices fixed
  (and corresponding spherical dual / I-minimal surface).

Each script is self-contained and can be executed independently.

## Requirements

- Julia ≥ 1.10  
- The required Julia packages are listed at the top of each script.

## Usage

From the repository root, run for example:

```bash
julia discrete_isothermic.jl
