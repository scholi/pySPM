from setuptools import setup, find_packages
import re

def description():
    with open('description.rst') as f:
        return f.read()

def find_version():
    with open("pySPM/__init__.py",'r') as fp:
        src = fp.read()
        version_match = re.search(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]", src, re.M)
        if version_match:
            return version_match.group(1)
        raise RuntimeError("Unable to find vesrion string.")
setup(
    name="pySPM",
    version=find_version(),
    description="library to handle SPM and ToF-SIMS data",
    long_description=description(),
    url="https://github.com/scholi/pySPM",
    author = "Olivier Scholder",
    author_email = "o.scholder@gmail.com",
    license="Apache 2.0",
    keywords='tof sims iontof spm sfm sxm afm kpfm pca imaging ita itm bruker nanonis',
    packages=find_packages(exclude=['contrib','docs','tests']),
    package_data={'pySPM':['data/elements.db','data/test.sxm']},
    include_package_data=True,
    install_requires=['numpy','scipy','pandas','scikit-image','scikit-learn','matplotlib'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)
