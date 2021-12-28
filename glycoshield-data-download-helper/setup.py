import os
import subprocess
from setuptools import setup
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install

def zenodo_downloader():
    doi="10.5281/zenodo.5337276"
    out=os.path.join(os.path.expanduser("~"), "GLYCAN_LIBRARY")
    if not os.path.exists(out):
        os.makedirs(out)
    cmd="zenodo_get --output-dir={} {} &>{}/zenodo_get.log".format(out, doi, out)
    subprocess.check_call(cmd, shell=True)

class develop(_develop):
    def run(self):
        zenodo_downloader()
        _develop.run(self)

class install(_install):
    def run(self):
        zenodo_downloader()
        _install.run(self)

setup(
    name='glycoshield-data-download-helper',
    version='0.1',
    packages=[],
    cmdclass={
        'develop': develop,
        'install': install,
    },
)

