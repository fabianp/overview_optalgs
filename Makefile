

s3_upload:
	aws s3 sync --region eu-west-1 --acl public-read public s3://fa.bianp.net/teaching/2018/eecs227at/
