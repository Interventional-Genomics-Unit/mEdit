# Installed Modules
import boto3


def main():
	s3_resource = boto3.resource('s3')

	def list_objects_in_bucket(bucket_name, prefix=''):
		# Create an S3 client
		s3 = boto3.client('s3')

		# List objects in the bucket with the specified prefix
		response = s3.list_objects_v2(Bucket=bucket_name, Prefix=prefix)

		# Iterate through the objects
		for obj in response.get('Contents', []):
			# Print the object key (file path)
			print(obj['Key'])

			# Check if the object is a directory (common S3 delimiter is '/')
			if obj['Key'].endswith('/'):
				# Recursively list objects in the subdirectory
				list_objects_in_bucket(bucket_name, prefix=obj['Key'])

	# Replace 'your-bucket-name' with your actual S3 bucket name
	list_objects_in_bucket(bucket)

	s3_resource.Object(first_bucket_name, first_file_name).download_file(
		f'/tmp/{first_file_name}')

if __name__ == "__main__":
	main()
