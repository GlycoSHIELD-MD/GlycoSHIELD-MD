from setuptools import setup
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install
from setuptools.command.test import test as _test

from download_helper import make_directory, download, unpack

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
        download()
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

