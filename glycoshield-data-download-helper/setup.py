import os
import sys
import glob
import shutil
import tempfile
import subprocess
from setuptools import setup
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install
from setuptools.command.test import test as _test


doi="10.5281/zenodo.5337276"
git="https://gitlab.mpcdf.mpg.de/MPIBP-Hummer/glycoshield-md-data.git"
out=os.path.join(os.path.expanduser("~"), "GLYCAN_LIBRARY")
tmp=os.path.expanduser("~")
pwd=os.getcwd()


def make_directory():
    if not os.path.exists(out):
        os.makedirs(out)

def zenodo_download():
    cmd="zenodo_get --output-dir={} {}".format(out, doi)
    subprocess.check_call(cmd, shell=True)

def gitlab_download():
    git_dir=tempfile.mkdtemp(dir=tmp)
    try:
        cmd="git clone --depth=1 {} {}".format(git, git_dir)
        subprocess.check_call(cmd, shell=True)
        os.path.join(git_dir, doi)
        file_list=glob.glob(os.path.join(os.path.join(git_dir, doi), "*"))
        for file_name in file_list:
            shutil.move(file_name, out)
    except:
        raise
    finally:
        shutil.rmtree(git_dir)

def download():
    gitlab_download()

def unpack():
    try:
        os.chdir(out)
        file_list=glob.glob("*.tar.gz")
        for file_name in file_list:
            cmd="tar -xf {}".format(file_name)
            subprocess.check_call(cmd, shell=True)
            os.remove(file_name)
    except:
        raise
    finally:
        os.chdir(pwd)

class develop(_develop):
    def run(self):
        make_directory()
        download()
        unpack()
        _develop.run(self)

class install(_install):
    def run(self):
        make_directory()
        download()
        unpack()
        _install.run(self)

class test(_test):
    def run(self):
        make_directory()
        #zenodo_download()
        gitlab_download()
        unpack()

setup(
    name='glycoshield-data-download-helper',
    version='0.1',
    packages=[],
    cmdclass={
        'develop': develop,
        'install': install,
        'test': test,
    },
)

