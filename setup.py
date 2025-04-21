# setup.py
from setuptools import setup, find_packages  # Import find_packages

ps = [
    "joblib", "numpy", "seaborn", "pandas", "scipy", "scikit-learn",
    "matplotlib", "tqdm", "pyBigWig", "click"
]

setup(
    name='ryder',
    version="1.0.0",
    author='Yaqiang Cao',
    author_email='caoyaqiang0410@gmail.com',
    url='https://github.com/YaqiangCao/ryder',
    keywords=
    'Epigenome, Cross-sample Normalization, ChIP-seq/DNase-seq/ATAC-seq/MNase-seq, Internal Reference',
    description=
    "Epigenome normalization with internal reference and variable features identifications.",
    classifiers=[
        'Environment :: Console',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    project_urls={
        'Source': 'https://github.com/YaqiangCao/ryder',
    },
    packages=find_packages(exclude=['test']),
    scripts=["src/paw.py", "src/patrol.py"],
    setup_requires=ps,
    install_requires=ps,
)
