from setuptools import setup

setup(
    name='GlycoSHIELD',
    url='http://glycoshield.eu/',
    version='0.1',
    packages=[
        'glycoshield',
    ],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'mdanalysis',
    ],
)

