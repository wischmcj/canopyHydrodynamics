# Getting started

## Installation

1. **Create a Virtual Environment**: elow wBe use the native 'venv' module to create a virtual environment. This is not strictly necessary, but it is a good practice to keep your project dependencies separate from your system dependencies in a virtual environment.
   ```bash
   python -m venv canHydroVenv
   ```

2. **Activate Your Environment**: Activate your virtual environment to install dependencies and run the project. The commands to activate the virtual environment depend on your operating system and shell. Below are the commands for activating the virtual environment in different operating systems and shells.

```bash
  # if using bash (Mac, Unix)
  source canHydroVenv/bin/activate
  # if using PowerShell (Windows)
  source canHydroVenv\Scripts\activate.ps1
```

3. **Install canoPyHydro**: canoPHydro is published with PyPA (the python packacing authority), so you can install the latest stable release of canoPyHydro using pip. This installs our latest stable release as well as several libraries required for the use of the package's features. canoPyHydro currently supports Python versions 3.9. 3.10 and 3.11.

```bash
   pip install canoPyHydro
```

You can also install also the latest development version by cloning the GitHub repository and using pip to install from the local directory:

```bash
$ pip install git+https://github.com/wischmcj/canopyHydrodynamics.git
```

4. **Set Configuration Options**: The default configuration file can be found at '/CanopyHydrodynamics/canopyhydro_config.toml'. Configuration options can be set by altering the contents of that file in place. Refer to the [configuration section below](https://canopyhydrodynamics.readthedocs.io/en/latest/getting_started.html#configuration) for more information on configuration options.

That's it! You're ready to start using canoPyHydro. You can find a basic tutorial below, as well as some more in depth use case examples in the [documentation](https://canopyhydrodynamics.readthedocs.io/en/latest/index.html).

# Configuration

For a quick start, the below can be added to the beginning of your script to automatically set the required configuration options.

```{python}
  from canopyhydro.configuration import *
```

However we recomend that you read on to understand the configuration options and how to set them.

## Environment Variables
The following environment variables are used to set the location of your configuration and logging files. These can be set or from in the terminal before running a script.

```bash
export CANOPYHYDRO_CONFIG=/path/to/your/config/file.toml
export CANOPYHYDRO_LOG_CONFIG=/path/to/your/logging/config/file.yml
```

or by placing the following code at the start of your script:

```{python}
  from canopyhydro.configuration import *
  config_file = os.environ["CANOPYHYDRO_CONFIG"] = f"{os.getcwd()}/canopyhydro_config.toml"
  log_config = os.environ["CANOPYHYDRO_LOG_CONFIG"] = f"{os.getcwd()}/logging_config.yml"
```

The default configuration file can be be found at '/CanopyHydrodynamics/canopyhydro_config.toml' and configuration options can be set by altering the contents of this file in place, or by setting the environment variables as described above to specify your own file location.

## QSM File structure
The [qsm] section details the column numbers in which each variable is stored in the input file. To read in this file correctly,
the following columns are required, and each row in the QSM must have a corresponding value for each:
        - cyl_id
        - parent_id
        - x (specify two columns as an array, one for x0 another for x1)
        - y (specify two columns as an array, one for y0 another for y1)
        - z (specify two columns as an array, one for z0 another for z1)
        - radius
        - volume
        - length
the below additional columns are required as well but can be left blank in the input file:
        - branch_order
        - reverse_branch_order
        - segment_id

## Model Parameters
More detail on these parameters can be found in the [documentation](). For now, all you need to ensure is that both
variables have an integer value
- min_len_drip_flow
    - must be > 0
- min_flow_grade_lim
    - must be > -pi/2 and

## Directories

The `root_dir` configuration variable is used to specify the root directory of the project. It is the base directory where all other directories are located.

The `input_dir` configuration variable is used to specify the directory where input files are located. This directory is where the library will look for QSM .csv files to read in.

The `output_dir` configuration variable is used to specify the directory where output files will be generated. This is the directory where the project will write results files (i.e. statistics, flow data) as well as pickle files and saved figures.

The `test_input_dir` configuration variable is used to specify the directory where test input files are located. This directory contains the files that are used for testing purposes, such as test data or sample input files. This directory is used by the test suite to run tests on the library, but will not be used in the main project.
