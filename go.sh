#!/bin/bash
#hugo
echo " "
#echo "=====Start==================="
echo "Push to github...."
echo " "
git add .

msg="Update the package `date`"
if [ $# -eq 1 ]
  then msg="$1"
fi
git commit -m "$msg"

# Push source and build repos.
git push origin master

# Come Back

echo "========== Done========"
echo " "

