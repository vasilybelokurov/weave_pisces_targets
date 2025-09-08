from __future__ import print_function
import os
from setuptools import setup
import glob
import subprocess


def get_revision():
    """
    Get the git revision of the code

    Returns:
    --------
    revision : string
        The string with the git revision
    """
    try:
        tmpout = subprocess.Popen('cd ' +
                                  os.path.dirname(os.path.abspath(__file__)) +
                                  ' ; git log -n 1 --pretty=format:%H',
                                  shell=True,
                                  bufsize=80,
                                  stdout=subprocess.PIPE).stdout
        revision = tmpout.read().decode()[:6]
        return revision
    except:
        return ''


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


VERSIONPIP = read('version.txt').rstrip()
VERSION = VERSIONPIP + '+dev' + get_revision()

with open('py/weave_chervin/_version.py', 'w') as fp:
    print('version="%s"' % (VERSION), file=fp)

setup(
    name="weave_chervin",
    version=VERSION,
    author="Sergey Koposov",
    author_email="skoposov _AT_ ed _DOT_ ac _DOT_ uk",
    description=("WEAVE selection"),
    license="BSD",
    keywords="WEAVE astronomy selection catalog survey",
    url="http://github.com/segasai/weave_galr",
    packages=[
        'weave_chervin', 'weave_chervin/selection', 'weave_chervin/utils'
    ],
    scripts=[fname for fname in glob.glob(os.path.join('bin', '*'))] +
    [fname for fname in glob.glob(os.path.join('scripts', '*sh'))],
    package_dir={'': 'py/'},
    package_data={
        'weave_chervin': [
            'data/WS2023B2-010_CatalogueTemplate.fits', 'conf/selection.yaml',
            'conf/target.yaml', 'conf/general.yaml', 'data/aps_flag.dat',
            'data/SharedSelectionMask.csv'
        ]
    },
    #    include_package_data=True,
    long_description=read('README.md'),
    install_requires=open('requirements.txt').readlines(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
