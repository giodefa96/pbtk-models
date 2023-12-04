# Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology

This project is based on the scientific paper "Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology" by Kannan Krishnan and Thomas Peyret. The paper provides quantitative descriptions of absorption, distribution, metabolism, and excretion of chemicals in biota.
Overview

## Overview  

PBTK models are increasingly being used as an effective tool for designing toxicology experiments and for conducting extrapolations essential for risk assessments. This project aims to provide a Python implementation of PBTK models for different species, including rats, humans, and now, trout.
Usage

To use this project, you can create a PBTK model for a specific species like this:

## Usage  

```python 
from pbtk import RatPBTKModel, TroutPBTKModel, Tetrachloroethane  
  
 
# For Rat  
rat = Rat()
molecule = Tetrachloroethane()
rat_model = RatPBTKModel(rat,molecule)    
rat_model.calculate_concentrations()    
rat_model.plot_results()    
rat_model.save_plots('rat_plot.png')    
```

```python  
# For Trout with Tetrachloroethane  
trout = Trout()  
tetrachloroethane = Tetrachloroethane()  
trout_model = TroutPBTKModel(trout, tetrachloroethane)  
trout_model.calculate_concentrations()    
trout_model.plot_results()    
trout_model.save_plots('trout_plot.png')    
```
 
This will create a PBTK model for a rat and a trout, calculate the concentrations, plot the results, and save the plots to files. If you want to add more species or molecules, you can create more subclasses of AbstractPBTKModel and Molecule respectively.
References

## References
 
Krishnan, K., & Peyret, T. (2009). Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology. In J. Devillers (Ed.), Ecotoxicology Modeling, Emerging Topics in Ecotoxicology: Principles, Approaches and Perspectives 2. Springer Science+Business Media, LLC.
Authors and acknowledgment

## Authors and acknowledgment 

Krishnan, K., Peyret, T. (2009). Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology. In: Devillers, J. (eds) Ecotoxicology Modeling. Emerging Topics in Ecotoxicology, vol 2. Springer, Boston, MA. 