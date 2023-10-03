#!/bin/bash 

# uncomment for running on a/the cluster
module load anaconda3/2022.05

# if ran with bash, conda init is necessary, if ran with source, conda init is not necessary
conda init bash

# create if command for when rerunning the script so it doesn't create the environment every time
if { conda env list | grep 'download_data'; } >/dev/null 2>&1; then
conda activate download_data
else
conda create -n download_data --yes -c conda-forge awscli
fi

# the while will loop through the string given to getops and
# for every letter/option it sees it will do something different
# this means the options can be given in any order.
# the options are stored in a variable named "flag"
# and OPTARG is the value of each variable at every run of the loop, ie each letter option
# the colons are used to indicate that each option is mandatory
# this is standard syntax for a case loop

while getopts ':i:k:r:f:u:d:' flag; do
	case "$flag" in
		i) 
			ACCESS_ID=${OPTARG};;
		k) 
			ACCESS_KEY=${OPTARG};;
		r)	
			REGION=${OPTARG};;
		f)
			OUTPUT=${OPTARG};;
		u)
			URL=${OPTARG};;
		d)
			OUT_DIR=${OPTARG};;
	esac
done

# create and edit the aws profile to be used
# a profile is simply a collection of aws access parameters, like username, password, file location etc
aws configure set aws_access_key_id $ACCESS_ID --profile download
aws configure set aws_secret_access_key $ACCESS_KEY --profile download
aws configure set region $REGION --profile download
aws configure set output $OUTPUT --profile download



aws s3 cp $URL $OUT_DIR --recursive --profile download

conda deactivate

# AKIA45O7KHZWD3U3PBF4
# FtIvAgOXYjyELRoqyo/BmZx57KUfxbsslFzxjd2U
# eu-north-1
# s3://bmk-az290.bmk-az741/bmk-az290/
# /share/pool/CompGenomVert/Bmk_RawData_Batch1/Owls/

# example run 1
# source ./download_data_aws_v2.sh -i "AKIA45O7KHZWD3U3PBF4" -k "FtIvAgOXYjyELRoqyo/BmZx57KUfxbsslFzxjd2U" -r "eu-north-1" -f "text" -u "s3://bmk-az290.bmk-az741/bmk-az290/" -d "/share/pool/CompGenomVert/Bmk_RawData_Batch1/Owls/"
# example run 2
# source ./download_data_aws_v2.sh -i "AKIA45O7KHZWD3U3PBF4" -k "FtIvAgOXYjyELRoqyo/BmZx57KUfxbsslFzxjd2U" -r "eu-north-1" -f "text" -u "s3://bmk-az290.bmk-az741/bmk-az290/Data/.check-ignore-list" -d "."

# fix line endings:
# sed -i 's/\r//' test1.sh
