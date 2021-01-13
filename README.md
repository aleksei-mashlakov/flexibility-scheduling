# flexibility-scheduling

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
