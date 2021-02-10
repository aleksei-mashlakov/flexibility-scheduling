# Prosumer flexibility scheduling

The codes and data provided here are used for the experiment simulations described in:

```
@article{mashlakov2021decentral,
  title={Decentralized cooperative scheduling of prosumer flexibility under forecast uncertainty},
  author={Mashlakov, A., Pournaras, E., Nardelli, P.H.J., and Honkapuro, S.},
  journal={},
  volume={},
  number={},
  pages={},
  year={2021},
  publisher={}
}
```

## Repository structure
<!--toc-->
 - [Data sets for forecasting and scheduling](#data)                    
 - [Prosumer flexibility modeling in R](#prosumer-flexibility-modeling-in-r)        
 - [Flexibility coordination package in java](#flexibility-coordination-package-in-java)    
 - [Battery control simulation in python](#battery-control-simulation)
<!--toc_end-->

#### Installation for Linux machine

    $ mkdir flexibility-scheduling
    $ cd flexibility-scheduling
    $ git clone https://github.com/aleksei-mashlakov/flexibility-scheduling.git .

## Data

**Carbon intensity** data is fetched from [National Grid ESO](https://data.nationalgrideso.com/carbon-intensity1/historic-generation-mix)

To get the data:

    $ wget https://data.nationalgrideso.com/backend/dataset/88313ae5-94e4-4ddc-a790-593554d8c6b9/resource/f93d1835-75bc-43e5-84ad-12472b180a98/download/df_fuel_ckan.csv -P ./datasets

The preprocessing was done using:

    $ python ./data_preprocessing/preprocess_carbon_intensity.py

**Household net load** data is fetched from [Network Revolution](http://www.networkrevolution.co.uk/resources/project-data/) project

To get the data:

    $ wget http://www.networkrevolution.co.uk/go.php?id=409&link=TC5.zip -P ./datasets
    $ wget http://www.networkrevolution.co.uk/go.php?id=409&link=TC2Auto.zip -P ./datasets
    $ unzip TC5.zip TC2Auto.zip -d ./datasets

The preprocessing was done using:

    $ python preprocess_net_load_data.py

**Grid Frequency** data is fetched from [National Grid ESO](https://data.nationalgrideso.com/system/system-frequency-data)

To get the data:

    $ bash ./datasets/get_frequency.sh

## Prosumer flexibility modeling in R

#### Prerequisites

**1. Installed R** (version>=3.6)

Instructions for linux are [here](https://blog.zenggyu.com/en/post/2018-01-29/installing-r-r-packages-e-g-tidyverse-and-rstudio-on-ubuntu-linux/) or [here](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-18-04-quickstart)

Ubuntu 18 in short:

    $ sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
    $ sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
    $ sudo apt-get update
    $ sudo apt-get install r-base

**2. (Optional) Installed Gurobi optimizer**

Gurobi speeds up the simulation but for the simple example we use the other solver.

Get and activate academic license. Check instructions from [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html).

    $ tar xvfz gurobi9.0.2_linux64.tar.gz
    $ sudo mv gurobi.lic /opt/gurobi902/linux64

Consider [this](https://stackoverflow.com/questions/44007425/gurobi-package-does-not-load-in-ubuntu-14-04-error-in-dyn-loadfile-dllpath) to run Gurobi in R.
Users of the bash shell should add the following lines to their .bashrc files in /etc/bash:

    export GUROBI_HOME="/opt/gurobi902/linux64"
    export PATH="${PATH}:${GUROBI_HOME}/bin"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

**3. Installed nohub**

    $ sudo apt-get install nohup

### Flexibility forecasting

Forecasts quantiles for carbon intensity, household net load, and battery grid frequency response.

    $ gunzip ./datasets/whole home power import 30 min 12 months averaged.csv.gz
    $ nohup Rscript ./R/plan_forecasting.R &


### Flexibility modeling

    $ nohup Rscript ./R/plan_modeling.R &


## Flexibility coordination package in java

#### Prerequisites

**1. Installed java**

    $ sudo apt install openjdk-11-jre-default

**2. Installed I-EPOS**

    $ wget https://github.com/epournaras/EPOS/releases/download/0.0.2/Release-0.0.2.zip -P ./I-EPOS/
    $ sudo apt install unzip
    $ unzip ./I-EPOS/Release-0.0.2.zip -d ./I-EPOS/

    <!--- to simply run jar: --->
        java -jar IEPOS-Tutorial.jar
        !the folder name should not contain gaps
    <!--- add Class-Path: . to MANIFEST file in IEPOS-Tutorial.jar (not necessary true) --->
        sudo chmod -x /prepare.sh
        sudo chmod -x /run_test.sh
    <!--- then use where 1 -- is an index of start dataset and 64 -- end dataset --->
        ./prepare.sh 1 64



## Battery control simulation
To run the experiment simulations:

    python ./Simulation/battery_control.py


























    <!--- To install R --->
    https://blog.zenggyu.com/en/post/2018-01-29/installing-r-r-packages-e-g-tidyverse-and-rstudio-on-ubuntu-linux/

    sudo apt-get install libeigen3-dev
    sudo apt-get install libmpfr-dev
    sudo apt-get install r-cran-rcppeigen

    <!--- To run Gurobi in linux --->
    tar xvfz gurobi9.0.2_linux64.tar.gz

    <!--- Users of the bash shell should add the following lines to their .bashrc files in /etc/bash: --->

    export GUROBI_HOME="/opt/gurobi902/linux64"
    export PATH="${PATH}:${GUROBI_HOME}/bin"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

    grbgetkey 68cfaeda-acf0-11ea-8319-0a7c4f30bdbe
    sudo mv gurobi.lic /opt/gurobi902/linux64
    # check gurobi_cl

    <!--- To run Gurobi in R:  https://stackoverflow.com/questions/44007425/gurobi-package-does-not-load-in-ubuntu-14-04-error-in-dyn-loadfile-dllpath --->
    cd /etc/ld.so.conf.d
    sudo nano x86_64-linux-gnu.conf
    ## add /opt/gurobi902/linux64/lib
    sudo ldconfig
    nohup Rscript ./script.R &

<!--- also consider --->

1.Add library path for R: Edit etc/R/ldpaths - open as sudo and add the following:

    : ${R_LD_LIBRARY_PATH=${GUROBI_HOME}/lib}

    #############################
    : ${R_LD_LIBRARY_PATH=${GUROBI_HOME}/lib}                -- this?
    : ${JAVA_HOME=/usr/lib/jvm/default-java}
    : ${R_JAVA_LD_LIBRARY_PATH=${JAVA_HOME}/lib/server}
    if test -n "/usr/lib/x86_64-linux-gnu"; then
    : ${R_LD_LIBRARY_PATH=${GUROBI_HOME}/lib}                -- this?


ALSO! Need to add a file in the directory etc/profile.d. Create file with the following text in it:

    export GUROBI_HOME="/opt/gurobi902/linux64"
    export PATH="${PATH}:${GUROBI_HOME}/bin"
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

and save as "gurobi902.sh" in that directory.
