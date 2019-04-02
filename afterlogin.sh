#!/bin/bash

apt-get install locales
locale-gen en_US.UTF-8
update-locale LANG=en_US.UTF-8
export LANG=en_US
export LC_ALL=en_US.UTF-8

cp -v /supp/*.py ./
cp -v /supp/_vimrc ./.vimrc
ln -sf /data/v2015

