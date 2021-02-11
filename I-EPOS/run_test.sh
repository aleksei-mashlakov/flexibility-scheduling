#!/bin/bash

# This shell runs IEPOS with following important parameters:
# - reads initial tree structure from datasets/initial-full-tree-<dataset>.csv file. Initial tree structure is given as already permuted list of nodes in this file.
#   However, which exact permutation is to be used is governed by parameter row >= 1.
# - after every tree reorganization, first plan that is selected (in iteration immediately following the reorganization) is presaved from some iteration
#   before reorganization. This 'some iteration' is given as offset parameter >= 1.
# - every experiment is repeated 100 times with different RANDOM SEED that will be used in random generator that randomly reorganizes the tree structure
# - iterations = 120, ON_CONVERGENCE strategy, 1000 agents, each with 10 plans, tree is binary

repetitions=1 #100
#numIterations=10
#max_offset=2 			# 1, 2, .., max_offset
#max_row=2 				  # 1, 2, .., max_row
#arguments="NOTHING"
max_beta=2
# beta_step=$(( 1 / $max_beta)) | sed 's/..$/.&/'
# echo "The beta step is "$beta_step
#datasets=("1") # "energy" "gaussian")

betas=(0.9995 1.0)
#betas=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.96 0.97 0.98 0.99 0.991 0.992 0.993 0.994 0.995 0.996 0.997 0.998 0.999 0.9991 0.9992 0.9993 0.9994 0.9995 0.9996 0.9997 0.9998 0.9999 0.99991 0.99992 0.99993 0.99994 0.99995 0.99996 0.99997 0.99998 0.99999 1.0)
# betas=(0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 \
#        0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 \
#        0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 \
#        0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 \
#        0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 \
#        0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 \
#        0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.7 \
#        0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.8 \
#        0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 \
#        0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.0)

# no params
run_jar() {
	java -jar "IEPOS-Tutorial.jar" # ${arguments}
}

run_all_repetitions() {
	for (( rep=0; rep<${repetitions} ; rep++ ))
	do
		#make_command $1 $2 #$3 ${randomSeed[$rep]}
		echo "DATASET: "$1", BETA: "$2 #", OFFSET: "$3", REP: "${rep}
		run_jar
	done
	dst_dir="output-"$1"-beta-"$2
	mv "output" $dst_dir
	rm -r log
}

run_all_betas() {
	for (( beta_idx=0; beta_idx<=${max_beta} ; beta_idx++ ))
	do
    beta=${betas[$beta_idx]}
    #alpha=${betas[$beta_idx]}
    echo $beta
    sed -i "45s/.*/weightsString = \'0.00,${beta}\'/" ./conf/epos.properties
		run_all_repetitions $1 $beta
	done
}

run_dataset() {
	sed -i "3s/.*/dataset=${1}/" ./conf/epos.properties
  run_all_betas $1
}

run_all_datasets() {
	for (( data_id=1; data_id<=${1} ; data_id++ ))
	do
		run_dataset $data_id
	done
}
echo "Dataset run "$1
#run_all_datasets $1
run_dataset $1
