from distutils.core import setup

setup(
    name="PyESPER",
    version="1.0.0",
    description="Python version of ESPERv1",
    author="LMD",
    author_email="lmdias@uw.edu",
    packages=["PyESPER"],
    install_requires=["numpy", "seawater", "scipy", "matplotlib", "PyCO2SYS"],
)
