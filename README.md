
This repository houses the python code relating to the below paper, and will recieve periodic updates as our research progresses.
A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
(https://www.researchgate.net/publication/375530854)

## Contents:

The interactive jupyter notebook under '.\Cylinders\cli.ipynb' displays the code written for the first draft of the above linked paper how it was run and reviewed. The remaning (majority) of this repository represents code written in the process of improving and productionalizing that code.

## Improvements since publication of preprint

The code that you see zipped here is representiitive of the ongoing work in our ***productionalizing*** branch, where you will see in progress improvements such as:

+ Functionality refactored into testable methods
+ Flexible, fast visualization functionality added
+ Linter(s) added for formatting and best practices adherence
+ Fully fledged logging functionality added for visibility into results
+ toml configuration added for model parameterization
+ A pytest based testing suite to surface potential disagreements with previous results 

## Development Setup:

- Pre-requisites

  1. Python version 3.9 or higher
  2. A Virtual environment
     - python -m venv ~\Venvs\CanopyHydroEnv
  3. Activate venv with 
      -  ~\Venvs\CanopyHydroEnv\Scripts\activate.ps1 (PowerShell)
  4. Install requirements 
      - pip install -r requirements_dev.txt
  5. Enabling Pre-commit lining
     2. Install pre-commit - pre-commit install

## Known Issues

- The projection algorithm is an approximation
  - e.g. a vector forming a 3,5,6 triangle (vector from (1,1,1)(4,6,7)) has an angle of 45 degrees or 0.785 rad with the XY plane but the algorithm returns .799 rad
- Cylinders occasionally have more than one possible drip point. Our algorithm chooses just one of theses and lists that as the drip node
  - Note that these possible drip nodes more often than not fall very close together
