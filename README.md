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

### Repository structure
<!--toc-->
    ├── [Data](#data)              # data sets for forecasting, scheduling             
    ├── [I-EPOS](#i-epos)          # flexibility coordination package in java
    ├── [R](#r)                    # prosumer flexibility modeling in R
    └── [Simulation](#simulation)  # article experiment simulation in python
<!--toc_end-->

### Data

    **Carbon intensity** data is retrieved from https://data.nationalgrideso.com/carbon-intensity1/historic-generation-mix#
    To get and pre-process the data:
    wget https://data.nationalgrideso.com/backend/dataset/88313ae5-94e4-4ddc-a790-593554d8c6b9/resource/f93d1835-75bc-43e5-84ad-12472b180a98/download/df_fuel_ckan.csv -P ./data
    python process_carbon_intensity.py

    **Household net load** data is retrived from http://www.networkrevolution.co.uk/resources/project-data/
    To get and pre-process the data:
    wget http://www.networkrevolution.co.uk/go.php?id=409&link=TC5.zip -P ./data
    wget http://www.networkrevolution.co.uk/go.php?id=409&link=TC2Auto.zip -P ./data

    **Grid Frequency** data is retrieved from https://data.nationalgrideso.com/system/system-frequency-data
    To get and pre-process the data:
    curl -L https://data.nationalgrideso.com/backend/dataset/cb1cc925-ecd8-4406-b021-3a3f368196e1/resource/f0933bdd-1b0e-4dd3-aa7f-5498df1ba5b9/download/f-2019-1.zip
    ...
    curl -L https://data.nationalgrideso.com/backend/dataset/cb1cc925-ecd8-4406-b021-3a3f368196e1/resource/f0933bdd-1b0e-4dd3-aa7f-5498df1ba5b9/download/f-2019-12.zip

### Simulation
To run the experiment simulations:

    python ./Simulation/battery_control.py

### Installation for Linux machine
    mkdir flexibility-scheduling
    cd flexibility-scheduling
    git clone https://github.com/aleksei-mashlakov/flexibility-scheduling.git .
    wget https://github.com/epournaras/EPOS/releases/download/0.0.2/Release-0.0.2.zip -P ./I-EPOS/
    sudo apt install unzip
    sudo apt install openjdk-11-jre-default

<!--- to simply run jar: --->
    java -jar IEPOS-Tutorial.jar
    !the folder name should not contain gaps
<!--- add Class-Path: . to MANIFEST file in IEPOS-Tutorial.jar (not necessary true) --->
    sudo chmod -x /prepare.sh
    sudo chmod -x /run_test.sh
<!--- then use where 1 -- is an index of start dataset and 64 -- end dataset --->
    ./prepare.sh 1 64

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
