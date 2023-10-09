# Native modules
import re
import os
# Installed Modules
import boto3


def list_objects_in_bucket(bucket_name, version, prefix='', continuation_token=None):
    objects_list = []

    # Create an S3 client
    s3 = boto3.client('s3')

    # Configure the request parameters, including the ContinuationToken
    params = {'Bucket': bucket_name, 'Prefix': prefix}
    if continuation_token:
        params['ContinuationToken'] = continuation_token

    # List objects in the bucket with the specified prefix
    response = s3.list_objects_v2(**params)

    # Iterate through the objects
    for obj in response.get('Contents', []):
        # Only compressed FASTA files are listed
        if re.search(r"{}.fa.gz$".format(version), obj['Key']):
            objects_list.append(str(obj['Key']))

    # Check if there are common prefixes (subdirectories)
    if 'CommonPrefixes' in response:
        for subdir in response['CommonPrefixes']:
            subdir_prefix = subdir['Prefix']
            # Recursively list objects in the subdirectory, passing the ContinuationToken
            subdir_objects, _ = list_objects_in_bucket(bucket_name, version, prefix=subdir_prefix, continuation_token=response.get('NextContinuationToken'))
            objects_list += subdir_objects

    return objects_list, response


def download_objects_from_bucket(bucket_name, objects_list, download_path):
	# Create an S3 client
    s3 = boto3.client('s3')
	for obj_key in objects_list:
		# Specify the local file path for downloading
		local_file_path = os.path.join(download_path, os.path.basename(obj_key))

		# Download the object from S3 to the local path
		s3.download_file(bucket_name, obj_key, local_file_path)
		# Optionally, you can print a message to indicate successful download
		print(f"Downloaded: {obj_key} to {local_file_path}")


def main():
	s3_bucket_name = 'human-pangenomics'
	s3_filepath = 'working/'
	assembly_version = "f1_assembly_v2_genbank"
	local_download_path = '/groups/clinical/projects/clinical_shared_data/hprc/aws/'
	# Replace 'your-bucket-name' with your actual S3 bucket name
	file_list, r = list_objects_in_bucket(s3_bucket_name, assembly_version, s3_filepath)

	# Download the objects to the specified path
	download_objects_from_bucket(s3_bucket_name, file_list, local_download_path)


if __name__ == "__main__":
	main()
