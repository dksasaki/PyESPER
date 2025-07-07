from distutils.core import setup

setup(
    name="mysci",
    version="1.0.0",
    description="Python ESPER",
    author="LMD",
    author_email="lmdias@uw.edu",
    packages=["PyESPER"],
    install_requires=["numpy", "seawater", "scipy", "matplotlib", "PyCO2SYS", "importlib"],
)
