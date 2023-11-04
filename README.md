# dripDropFlow
Houses code relating to 'Whence Stemflow &amp; Throughfall Drip Points?'
Reflects the authors' naiive view of the project at its outset

<h1>CanopyHydrodymaics</h1>

<h2>Setup:</h2>
- Pre-requisites
  1.  Python version 3.9 or higher
  2.  A Virtual environment
      - python -m venv ~\Venvs\CanopyHydroEnv

- Start-up
  1. Activate venv with -  ~\Venvs\CanopyHydroEnv\Scripts\activate
  2. Install requirements - pip install requirements.txt
  3. Install pre-commit - pre-commit install

- pre-commit
  1. set-up - https://pre-commit.com/
  2. run on all files - pre-commit run --all-files

- pre-commit
  1. set-up - https://pre-commit.com/
  2. run on all files - pre-commit run --all-files


<h2>Contents:</h2>

- "_src" Directories
  1. **python_src:** Refactors the code originally developed in matlab in python. Allows for the use of more robust (not to mention, free) libraries for network analysis (networkx)
  2. **matlab_src:** Original code, as of 4/27 has the additional benefit of also plotting the cyliders

- Data
  1. **:** contains sample files to allow for running of the program.
  2. **Cyl Fitting Types:**  EXTRAPOLATETODTM, ALLOMETRICGROWTHVOLUME, SPHEREFOLLOWING, MEDIAN

<h2>Commands:</h2>

- Activate venv
  - in Powershell
  - cd C:\Users\wisch\Documents\GitProjects\cylinders_personal
  - venv\scripts\activate.ps1
- Run the projects
  - load graphs from last files
    - python -m cylinders load
  - create graphs from cyl data
    - python -m cylinders read

<h2>Approach</h2>

- Pixel data from Lidar is read in to an object called cylinder
    - An array variable of said cylider (pSV) is the created that stores the projected vetors for each cylinder
    - The proejcted vecors (in conjunctuon with the parent-child relationships from the raw data) are then used to create an adjacency matrics
    - The adjacency matrix is then converted into a graph with networkx
    - This graph is then used to create an array? of stempaths for analysis
