# colibre-skirt
This repository contains a tutorial (`tutorial.ipynb`) for postprocessing a single galaxy from the [COLIBRE](https://colibre.strw.leidenuniv.nl/) simulations with the [SKIRT](https://skirt.ugent.be/root/_home.html) 3D dust radiative transfer code.
## Access to COLIBRE data
To access the COLIBRE data, you will need an account on the cosma8 machine. The tutorial notebook illustrates how to interface with the COLIBRE data for a simple SKIRT workflow and how to derive the necessary SKIRT input files. These input files are also included in this repository (`SKIRTinputFiles/`).
## SKIRT setup
The SKIRT setup corresponds mostly to the one described in Gebek et al. 2026. The configuration file is included in this repository (`template_v5.0.ski`). To install SKIRT, follow the instructions from the SKIRT [website](https://skirt.ugent.be/root/_home.html).
## Postprocessing more galaxies
This tutorial covers the basic COLIBRE-SKIRT worfklow for a single galaxy. To postprocess more galaxies, an efficient workflow will require iterating over galaxies and using HPC systems. The `hpc_workflow_anna` branch of this repository contains such an example workflow.
## Contact
In case of questions or comments, please reach out to Andrea Gebek (andrea.gebek@ugent.be).