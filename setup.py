from setuptools import setup, find_packages

setup(
    name='MsPAC',
    description='',
    packages=find_packages(),
    include_package_data=True,
    entry_points = {
        'console_scripts': ['MsPAC = MsPAC.MsPAC:main'],
        },
    platforms='any'
)
