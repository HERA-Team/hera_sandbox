#! /bin/bash
#update qmaster
sudo python setup.py install
#update the stills
sudo python setup.py install --prefix=/srv/precise_root.x86_64/usr/local/
