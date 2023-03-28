
from setuptools import setup, find_packages

requirements = open('requirements.txt').read().splitlines()

setup(
    name='dms_ci',
    version='1.0.0',
    url='https://github.com/rouskinlab/DMS_CI',
    author='Yves Martin',
    author_email='yves@martin.yt',
    description='A method to calculate confidence intervals for DMS-MaPseq data',
    packages=find_packages(),    
    install_requires=requirements,
)