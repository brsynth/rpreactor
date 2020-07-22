from setuptools import setup
import versioneer

requirements = [
    'rdkit==2019.03.1.0',
    'rpchemtools==0.4.0'
]

setup(
    name='rpreactor',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Apply reaction rules and parse results',
    license='MIT',
    author='Thomas Duigou',
    author_email='thomas.duigou@inrae.fr',
    url='https://github.com/tduigou/rpreactor',
    packages=['rpreactor'],
    install_requires=requirements,
    keywords='rpchemtools',
    classifiers=[
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
