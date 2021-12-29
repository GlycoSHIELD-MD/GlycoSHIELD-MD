#!/usr/bin/env python3


import os
import sys
import glob
import shutil
import tempfile
import subprocess


# --- configuration options ---
# DOI of the dataset as published on Zenodo
doi="10.5281/zenodo.5337276"
# Data repository on MPCDF GitLab
git="https://gitlab.mpcdf.mpg.de/MPIBP-Hummer/glycoshield-md-data.git"
# data folder for the glycan library, default is relative to the base directory of the 'glycoshield' package
out=os.path.join(
    os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")),
    "GLYCAN_LIBRARY")
# temporary directory prefix for git clone
tmp=os.path.expanduser("~")
# save initial pwd
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
        dest_list=[os.path.join(out, os.path.basename(file_name)) for file_name in file_list]
        for file_in, file_to in zip(file_list, dest_list):
            shutil.move(file_in, file_to)
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


if __name__ == "__main__":
    make_directory()
    download()
    unpack()

