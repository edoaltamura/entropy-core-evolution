# Entropy Evolution: Reproducing the Emergence of Entropy Plateaus

## Overview  
This repository contains all data and analysis scripts used to reproduce the 
results of Altamura et al. (2025, Paper II), which investigates the emergence of entropy 
plateaus in simulated galaxy groups and clusters using cosmological hydrodynamic 
zoom-in simulations with the SWIFT-EAGLE model. It builds directly on the $z = 0$ 
entropy profile analysis presented in Altamura et al. (2023, Paper I).

<details>
  <summary>Summary of main science</summary>

- Entropy plateaus emerge at characteristic halo-mass scales. Simulations of a galaxy group ($M_
  {500}\simeq8.8\times10^{12}\,M_\odot$) and a cluster ($M_{500}\simeq2.9\times10^{14}\,M_\odot$)
  show that once a halo reaches $M\sim10^{12}\,M_\odot$, its entropy profile flattens at the 
  virial radius. As the halo grows to $\sim10^{13}\,M_\odot$, the plateau extends inward, and by 
  $\sim10^{14}\,M_\odot$ a fully isentropic core is established.

- AGN feedback is the principal mechanism.
Lagrangian tracking of gas particles reveals that AGN outbursts expel low-entropy gas before it can accrete into the core, replacing it with higher-entropy material and erasing the central gradient needed for a cool core.

- Transition coincides with peak SMBH activity.
The onset of the entropy plateau at $M\sim10^{12}\,M_\odot$ aligns with the maximum in the 
  specific black-hole accretion rate, indicating a shift from supernova-dominated to AGN-dominated thermodynamic regulation.

- Numerical convergence.
High-resolution runs (gas particle mass $m_{\rm gas}\lesssim2.3\times10^5\,M_\odot$) confirm 
  that the entropy plateau persists even when subgrid physics is resolved on smaller scales.

- Comparison with observations.
XMM–*Newton* studies of local groups report entropy excesses and flat cores consistent with the 
  predicted plateaus, while many clusters still exhibit steep, cool-core power laws. Reproducing the observed diversity of entropy profiles remains a challenge.

- Implications for AGN subgrid modeling.
The tendency to over-eject low-entropy gas suggests that current feedback prescriptions may be too aggressive at group scales. Adaptive efficiency schemes or hybrid thermal–kinetic models may be required to recover the full spectrum of entropy shapes without compromising other cluster properties.
</details>

## Repository structure  
```text
├── data/               # Simulation data products
├── analysis_scripts/   # Analysis scripts to generate the data products
└── figures_scripts/    # Scripts to generate the figures from data products
```

## Installation  
1. Clone this repository
```bash
git clone https://github.com/edoaltamura/entropy-core-evolution.git
cd entropy-core-evolution
```

2. Create and activate the environment


## Citation
If you use this repository or its data in your work, please cite:

```text
```

## License
This project is licensed under the Apache License Version 2.0.