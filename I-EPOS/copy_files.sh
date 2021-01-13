dataset_ids_start=$1
dataset_ids=$2

get_current_path() {
	read base_path <<< $(pwd)
	base_path=$base_path"/"
	echo $base_path
}

make_and_copy_directories() {
	for (( dataset_id=$dataset_ids_start; dataset_id<=$dataset_ids ; dataset_id++ ))
	do
		dir_name="run_beta_day"$dataset_id
    #echo $dir_name
    src_dir=$base_path$dir_name"/output-"$dataset_id"-beta-0.9998"
    #echo $src_dir
    dst_dir="/home/ubuntu/puhti_scratch/flex_run/"$dir_name
    #echo $dst_dir
    mkdir $dst_dir
		cp -r $src_dir $dst_dir"/"
	done
}

get_current_path
make_and_copy_directories
