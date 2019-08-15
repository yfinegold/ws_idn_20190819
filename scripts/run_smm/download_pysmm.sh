#!/bin/bash
# Install PYSMM

echo "Installing PYSMM"

mkdir smm
cd smm
virtualenv -p python2.7 env --system-site-packages
source env/bin/activate
pip install pysmm
pip install ipyleaflet

wget -O pysmm_upgrade.zip https://goo.gl/JZLCo9
unzip -o pysmm_upgrade.zip -d env/lib/python2.7/site-packages/ 
python -m ipykernel install --user --name=env  
wget -O run_pysmm.ipynb https://goo.gl/NbpUWr
deactivate

echo "PYSMM installation complete"