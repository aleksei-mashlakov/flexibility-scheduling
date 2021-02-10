
## Prosumer flexibility modeling

**R scripts** contain:

- predict_garch.R: GARCH mean and covariance matrix forecasts.
- predict_quantiles.R: Generates temporal scenarios and derives marginal quantiles per time stamp.
- plan modeling.R: Multi-objective prosumer optimization

### Prerequisites

**1. Installed R** (version>=3.6)
Instructions for linux are [here](https://blog.zenggyu.com/en/post/2018-01-29/installing-r-r-packages-e-g-tidyverse-and-rstudio-on-ubuntu-linux/) or [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-18-04-quickstart)

Ubuntu 18:

    $ sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
    $ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    $ sudo apt-get update
    $ sudo apt-get install r-base


    $ sudo apt-get install libeigen3-dev
    $ sudo apt-get install libmpfr-dev
    $ sudo apt-get install r-cran-rcppeigen

**2. (Optional) Installed Gurobi optimizer**

Gurobi speeds up the simulation but for the simple example we use the other solver.

Get and activate academic license from [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html).

    $ tar xvfz gurobi9.0.2_linux64.tar.gz
    $ sudo mv gurobi.lic /opt/gurobi902/linux64

Consider [this](https://stackoverflow.com/questions/44007425/gurobi-package-does-not-load-in-ubuntu-14-04-error-in-dyn-loadfile-dllpath) to run Gurobi in R.

E.G., users of the bash shell should add the following lines to their .bashrc files in /etc/bash:

    export GUROBI_HOME="/opt/gurobi902/linux64"
    export PATH="${PATH}:${GUROBI_HOME}/bin"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

**3. Installed nohub**

    $ sudo apt-get install nohup

### Flexibility forecasting

Forecasts quantiles for carbon intensity, household net load, and battery grid frequency response

    $ nohup Rscript ./R/plan_forecasting.R &


    ### Flexibility forecasting

    $ nohup Rscript ./R/plan_modeling.R &
