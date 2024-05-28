"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
https://docs.python.org/3/distutils/setupscript.html#distutils-installing-scripts
"""
from setuptools import setup, find_packages


EXCLUDE_FROM_PACKAGES = [
    'vibronic_models',
    'tests',
]


def setup_package():

    setup_info = dict(
        name='vmio',
        version='1.0',
        python_requires='~={}.{}'.format(3, 9),
        url='https://github.com/ngraymon/vmio',
        author='Neil Raymond',
        author_email='neil.raymond@uwaterloo.ca',
        description='Assorted scripts and code for handling vibronic models and interfacing with MCTDH',
        packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
        include_package_data=True,
        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Chemistry',
        ],
        keywords='path_integral_monte_carlo quantum_mechanics chemistry vibronic',
        # https://packaging.python.org/en/latest/requirements.html
        install_requires=[
            'numpy>=1.12',
            'parse==1.8.2',
        ],
        extras_require={
            'dev': ['check-manifest'],
            'test': ['pytest', ],
        },
    )

    setup(**setup_info)


if __name__ == '__main__':
    setup_package()
