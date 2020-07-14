#!/bin/sh

if [ -z "$1" ]
  then
    echo "No check pattern argument supplied" >&2
    exit 1
fi

if [ -z "$2" ]
  then
    echo "No verify pattern argument supplied" >&2
    exit 1
fi

checkfiles=$1
infiles=$2
#echo $checkfiles
#echo $infiles

echo '\nCalculating check file hashes...'
md5sum $checkfiles

echo '\nCalculating input file hashes...'
md5sum $infiles

echo '\nComparing hash of file of hashes...'
checkver=$(md5sum $checkfiles | awk '{print $1}' | md5sum | awk '{print $1}')
echo $checkver

inver=$(md5sum $infiles | awk '{print $1}' | md5sum | awk '{print $1}')
echo $inver

if [ "$checkver" == "$inver" ]
then
  echo "Hashes match"
  exit 0
else
  echo "Hashes do not match" >&2
  exit 1
fi