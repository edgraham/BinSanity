from setuptools import setup
from codecs import open
from os import path
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
setup(
    name='Binsanity',    # This is the name of your PyPI-package.
    version='0.5.4',
    description='Method to cluster contigs based a biphasic method with coverage and composition',
    url="https://github.com/edgraham/BinSanity",
    author="Elaina Graham",
    author_email="egraham147@gmail.com",
    license="GPL3",
    install_requires=['biopython','pandas>=0.13.0','scipy>=0.13.0','checkm-genome','scikit-learn>=0.23'],

    scripts=['utils/identifyHMM','utils/concat','utils/simplify-fasta','utils/get-ids','bin/Binsanity-lc','bin/Binsanity','bin/Binsanity-refine','bin/Binsanity-wf',"bin/Binsanity-profile","utils/transform-coverage-profile","utils/bin_evaluation","utils/checkm_analysis","Binsanity2-Beta/Binsanity2-beta"]      ,
)

