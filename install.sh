# this was tested in Ubuntu 22.04
# I don't know if it will work in windows computers

conda create -n esper2 --yes
conda activate esper2

# do not install jupyter using conda install - jupyer won't work
# because pip will break its installation
# you probably can include cartopy here
conda install -c conda-forge numpy==1.16.0 xarray netCDF4 

git clone https://github.com/dksasaki/PyESPER.git
cd PyESPER

pip install --upgrade-strategy only-if-needed .
pip install jupyter
