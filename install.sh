#!/bin/bash
# script for installation of MGERT and test_dataset.tgz
# requires sudo

# check permissions
if [ "$EUID" -ne 0 ]
  then echo "Please run as root"
  exit
fi

mkdir /usr/share/mgert
chmod +x MGERT.py
cp MGERT.py /usr/local/bin
cp test_dataset.tgz /usr/share/mgert
printf "Installation completed.\nMGERT.py is installed in /usr/local/bin\ntest datasetis copied in /usr/share/mgert"
