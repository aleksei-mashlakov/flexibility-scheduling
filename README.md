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

To get the carbon intensity data:

    $ wget https://data.nationalgrideso.com/backend/dataset/88313ae5-94e4-4ddc-a790-593554d8c6b9/resource/f93d1835-75bc-43e5-84ad-12472b180a98/download/df_fuel_ckan.csv -P ./datasets

The preprocessing was done using:

    $ python ./data_preprocessing/preprocess_carbon_intensity.py

**Household net load** data is fetched from [Network Revolution](http://www.networkrevolution.co.uk/resources/project-data/) project

To get the net load data:

    $ wget http://www.networkrevolution.co.uk/go.php?id=409&link=TC5.zip -P ./datasets
    $ wget http://www.networkrevolution.co.uk/go.php?id=409&link=TC2Auto.zip -P ./datasets
    $ unzip TC5.zip TC2Auto.zip -d ./datasets

The preprocessing was done using:

    $ python preprocess_net_load_data.py

**Grid Frequency** data is fetched from [National Grid ESO](https://data.nationalgrideso.com/system/system-frequency-data)

To get the frequency data for year 2019:

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

    $ export GUROBI_HOME="/opt/gurobi902/linux64"
    $ export PATH="${PATH}:${GUROBI_HOME}/bin"
    $ export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

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

**2. Installed I-EPOS** (it is already in the repository but the full instructions below)

To run 1 day example of flexibility coordination. For the rest days you need to run **Flexibility modeling** and place the results in ./I-EPOS/datasets.

    $ cd ./I-EPOS
    $ ./prepare.sh 1 1

Full I-EPOS instructions below:

    $ wget https://github.com/epournaras/EPOS/releases/download/0.0.2/Release-0.0.2.zip -P ./I-EPOS/
    $ sudo apt install unzip
    $ unzip ./I-EPOS/Release-0.0.2.zip -d ./I-EPOS/
    $ rm -rf ./I-EPOS/Release-0.0.2.zip
    $ unzip ./I-EPOS/Release-0.0.2.zip -d ./I-EPOS/ && mv ./I-EPOS/Release-0.0.2/* ./I-EPOS/ && rm -rf ./I-EPOS/Release-0.0.2*
    $ mv ./datasets/1/ ./I-EPOS/datasets/1
    $ change number of agent to 150 and  ./I-EPOS/conf  
    $ change ./I-EPOS/conf/epos.properties file
    ...............................
    ### Dataset ###
    #The folder name in the datasets path. Make sure it has no spaces, tabs or newlines (alphanum and underscore preferred)
    dataset=1


    ### Basic epos properties ###
    # any integer > 0
    numSimulations=50

    # any integer > 0
    numIterations=30

    # any integer > 0
    numAgents=150

    # any integer > 0
    numPlans=19

    # any integer > 0
    numChildren=2

    # exact dimensionality from the dataset
    planDim=48
    ..................................

## Battery control simulation

The battery control simulations were conducted with the following file:

    $ python ./Simulation/battery_control.py
