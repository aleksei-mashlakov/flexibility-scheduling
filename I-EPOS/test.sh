batch_size=8
echo "Batch size is "$batch_size
dataset_ids_start=$1
dataset_ids=$2
echo "Dataset numbers are "$dataset_ids
num_batches=$(( (dataset_ids-dataset_ids_start-1) / batch_size))
echo "Number of batches is "$num_batches
batch_remainder=$((dataset_ids % batch_size))
echo "Batch remainder is "$batch_remainder

run_batches(){
  if (($batch_remainder > 0 )); then
    num_batches=$(($num_batches + 1))
    echo "Number of batches is "$num_batches
    #
  fi
  for (( b=1; b<=$num_batches ; b++ ))
  do
    end=$((dataset_ids_start-1 + b * batch_size ))
    start=$((dataset_ids_start-1 + (b-1) * batch_size + 1 ))
    echo $end
    echo $start
    #run_all_batches $start $end
  done
}

run_all_batches() {
	pids=()
	for (( dataset_id=$1; dataset_id<=$2 ; dataset_id++ ))
  do
    #echo $1
    echo $dataset_id

    if (( $dataset_id > $dataset_ids )); then
      echo "Max number of datasets is reached "$dataset_ids
      break;
    fi
		#dir_name="run_beta_day"$dataset_id
		#cd $dir_name
		#./${run_filename} ${dataset_id} ${jar_filename} &
		#echo $!
		#read pids <<< $!
		#cd ..
	done

	for (( i=$1; i<=$2 ; i++ ))
	do
		wait ${pids[$i]}
	done
}

run_batches
