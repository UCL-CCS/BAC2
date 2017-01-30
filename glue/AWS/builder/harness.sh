#! /bin/bash

source /home/ubuntu/.profile

/home/ubuntu/runbac.py > /home/ubuntu/builder.log 2>&1

OUT=$?
if [ $OUT -eq 0 ];then
   echo "DROPPING TO SHELL"
else
  bucket=$(/home/ubuntu/getbucket.pl)

  aws s3 cp /home/ubuntu/builder.log s3://$bucket/builder.log

  echo "FINISHED RUNNING, CLOSING INSTANCE"
  sudo shutdown -h now
fi
