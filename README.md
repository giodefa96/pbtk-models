# Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology  
  
This project is based on the scientific paper "Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology" by Kannan Krishnan and Thomas Peyret. The paper provides quantitative descriptions of absorption, distribution, metabolism, and excretion of chemicals in biota.   
  
## Overview  
  
PBTK models are increasingly being used as an effective tool for designing toxicology experiments and for conducting extrapolations essential for risk assessments. This project aims to provide a Python implementation of PBTK models for different species, including rats and humans.  
  
## Usage  
  
To use this project, you can create a PBTK model for a specific species like this:  
  
```python  
export PYTHONPATH="${PYTHONPATH}:src"
model = RatPBTKModel()  
model.calculate_concentrations()  
model.plot_results()  
model.save_plots('plot1.png', 'plot2.png')  
```

This will create a PBTK model for a rat, calculate the concentrations, plot the results, and save the plots to files. If you want to add more species, you can create more subclasses of AbstractPBTKModel.

## References

 
Krishnan, K., & Peyret, T. (2009). Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology. In J. Devillers (Ed.), Ecotoxicology Modeling, Emerging Topics in Ecotoxicology: Principles, Approaches and Perspectives 2. Springer Science+Business Media, LLC.

## Authors and acknowledgment


Krishnan, K., Peyret, T. (2009). Physiologically Based Toxicokinetic (PBTK) Modeling in Ecotoxicology. In: Devillers, J. (eds) Ecotoxicology Modeling. Emerging Topics in Ecotoxicology, vol 2. Springer, Boston, MA. https://doi.org/10.1007/978-1-4419-0197-2_6