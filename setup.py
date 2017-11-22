from setuptools import setup, find_packages

setup(
    name="pySPM",
    version='0.2',
    description="library to handle SPM and ToF-SIMS data",
    long_description="library to handle SPM and ToF-SIMS data",
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
        'Development Status :: 3 - Alpha',
        'Intended Audiance :: Developers',
        'Topic :: Scientific Data Processing :: Image processing',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)
