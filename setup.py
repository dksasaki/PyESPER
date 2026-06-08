from distutils.core import setup
import glob
import os

# Get all .py files in NeuralNetworks
neural_net_files = glob.glob("NeuralNetworks/*.py")

setup(
    name="PyESPER",
    version="1.0.0",
    description="Python version of ESPERv1",
    author="LMD",
    author_email="lmdias@uw.edu",
    packages=["PyESPER"],
    package_data={
        'PyESPER': ['NeuralNetworks/*.py', 'NeuralNetworks/__init__.py']
    },
    install_requires=["numpy", "seawater", "scipy", "matplotlib", "PyCO2SYS"],
)