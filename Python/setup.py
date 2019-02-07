from setuptools import setup

setup(
    name='airglowrsss',
    packages=['airglowrsss'],
    author="Brian Harding",
    description="ICON Project Files",
    license="GPLv3",
    url="https://github.com/bharding512/airglowrsss",
    install_requires=[
        "ephem",
        "matplotlib",
        "python-dateutil",
        "numpy",
        "netCDF4",
        "scipy",
        "pytz",
        "pandas",
        "Pillow",
    ],
    classifiers=[
        "Programming Language :: Python :: 2",
    ]
)
