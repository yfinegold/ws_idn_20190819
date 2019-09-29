#!/bin/bash
# Install PYSMM

echo "Installing PYSMM"
cd 
mkdir smm
cd smm
virtualenv -p python2.7 env --system-site-packages
source env/bin/activate
pip install pysmm
pip install ipyleaflet

wget -O pysmm_upgrade_20190922.zip https://www.dropbox.com/s/nwil00kw407ock4/pysmm_upgrade_20190922.zip
unzip -o pysmm_upgrade_20190922.zip -d env/lib/python2.7/site-packages/ 
python -m ipykernel install --user --name=env  
wget -O run_pysmm.ipynb https://goo.gl/NbpUWr
deactivate

echo "PYSMM installation complete"
