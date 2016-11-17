from setuptools import setup
from codecs import open
from os import path
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
setup(
    name='Binsanity',    # This is the name of your PyPI-package.
    version='0.1.8.3',
    description='Method to cluster contigs based a biphasic method with coverage and composition',
    url="https://github.com/edgraham/BinSanity",
    author="Elaina Graham",
    author_email="egraham147@gmail.com",
    license="GPL3",
    install_requires=['biopython','pandas>=0.13.0','scipy>=0.13.0','checkm-genome'],
    scripts=['Bin/Binsanity','Bin/Binsanity-refine','Bin/Binsanity-wf',"Bin/Binsanity-profile","Utils/transform-coverage-profile","Utils/bin_evaluation","Utils/checkm_analysis"]      ,            
)
