from setuptools import setup, find_packages
import versioneer

# WARNING: installation with setuptool or pip is not recommended since it does not list all dependencies.
# For instance, rdkit MUST be installed but WILL NOT be installed via pip. Prefer using conda.
# Have a look at the conda recipe (recipe/meta.yml) for a complete list of dependencies.

setup(
    name='rpreactor',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Handle biochemical reaction rules.',
    license='MIT',
    author='Thomas Duigou',
    author_email='thomas.duigou@inrae.fr',
    url='https://github.com/brsynth/rpreactor',
    packages=find_packages(),
    keywords=['rpreactor'],
    classifiers=[
        'Topic :: Scientific/Engineering',
    ]
)
