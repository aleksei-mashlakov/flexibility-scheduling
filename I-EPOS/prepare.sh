#!/bin/bash

# This script expects the following files and folders in the same directory (on the same level as this script):
#  - runnable jar IEPOS.jar
#  - runnable script run_test.sh
#  - folder 'conf' containing configuration files
#  - folder 'datasets' containing folder with plans. Three different datasets are supported:
#     - 'gaussian' -> requires command line argument 0
#     - 'bicycle'  -> requires command line argument 1
#	    - 'energy'   -> requires command line argument 2

# use batches when simulating for large number of days
batch_size=8 #$((nproc --all))
echo "Batch size is "$batch_size
dataset_ids_start=$1
dataset_ids=$2
echo "Dataset numbers are "$dataset_ids
num_batches=$(( (dataset_ids-dataset_ids_start-1) / batch_size))
echo "Number of batches is "$num_batches
batch_remainder=$(( (dataset_ids -dataset_ids_start) % batch_size))
echo "Batch remainder is "$batch_remainder

jar_filename="IEPOS-Tutorial.jar"
run_filename="run_test.sh"
conf_filename="conf"
base_path="NOTHING"

make_and_fill_directories() {
	for (( dataset_id=$dataset_ids_start; dataset_id<=$dataset_ids ; dataset_id++ ))
	do
		dir_name="run_beta_day"$dataset_id
		mkdir $dir_name
		cp -r $conf_filename $jar_filename $run_filename $dir_name
	done
}

get_current_path() {
	read base_path <<< $(pwd)
	base_path=$base_path"/"
	echo $base_path
}

make_dataset_symbolic_links() {
	for (( dataset_id=$dataset_ids_start; dataset_id<=$dataset_ids ; dataset_id++ ))
	do
		dir_name="run_beta_day"$dataset_id
    datasets_filename=$base_path"datasets"
    #datasets_filename="/home/ubuntu/R/plans_sonnen/" #$base_path"datasets"
		destination_path=$base_path$dir_name

		ln -s $datasets_filename $destination_path
	done
}

run_all() {
	pids=()
	for (( dataset_id=$dataset_ids; dataset_id<=$dataset_ids ; dataset_id++ ))
  do
		dir_name="run_beta_day"$dataset_id

		cd $dir_name

		./${run_filename} ${dataset_id} ${jar_filename} &
		echo $!
		read pids <<< $!
		cd ..
	done

	for (( i=$dataset_ids; i<=$dataset_ids ; i++ ))
	do
		wait ${pids[$i]}
	done
}

run_all_batches() {
	pids=()
	for (( dataset_id=$1; dataset_id<=$2 ; dataset_id++ ))
  do
    #echo $1
    #echo $dataset_id
    if (( $dataset_id <= $dataset_ids )); then
      #echo "Max number of datasets is reached "$dataset_ids
      #break;
  		dir_name="run_beta_day"$dataset_id
  		cd $dir_name
  		./${run_filename} ${dataset_id} ${jar_filename} &
  		echo $!
  		read pids <<< $!
  		cd ..
    fi
	done

	for (( i=$1; i<=$2 ; i++ ))
	do
    if (( $i <= $dataset_ids )); then
      wait ${pids[$i]}
    fi
	done
}

run_batches(){
  if (($batch_remainder > 0 )); then
    num_batches=$(($num_batches + 1))
    echo "Number of batches is "$num_batches
  fi

  for (( b=1; b<=$num_batches ; b++ ))
  do
    end=$(( dataset_ids_start-1 + b * batch_size ))
    start=$(( dataset_ids_start-1 + (b-1) * batch_size + 1 ))
    #echo $end
    #echo $start
    run_all_batches $start $end
  done
}

get_current_path
make_and_fill_directories
make_dataset_symbolic_links
#run_batches
run_all
echo "Done&Done."
