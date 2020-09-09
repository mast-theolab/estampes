import setuptools

with open('README.md', 'r') as fobj:
    long_description = fobj.read()
short_desc = 'A simple API to parse and process data for computational' +\
             ' spectroscopy'

setuptools.setup(
    name='estampes',
    version='0.1.2',
    author='Julien Bloino',
    author_email='julien.bloino@gmail.com',
    description=short_desc,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/jbloino/estampes',
    license='MIT',
    packages=setuptools.find_packages(include=['estampes', 'estampes.*',
                                               'progs.esparser']),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Visualization'
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'esparser=progs.esparser:main',
        ],
    },
)
