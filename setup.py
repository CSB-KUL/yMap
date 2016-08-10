from setuptools import setup, find_packages

setup(
	name = "ymap",
	version = "1.0",
	description = "An automated method to map yeast variants to protein modifications and functional regions",
	url = "https://github.com/CSB-KUL/yMap",
	author = "Ahmed Arslan and Vera van Noort",
	author_email = "vera.vannoort@kuleuven.be",
	maintainer = "Ahmed Arslan",
	maintainer_email = "ahmed.arslan@kuleuven.be",
	license = "MIT",
	classifiers=[ 
	'Intended Audience :: User ((Computational)Biologists)',
	'Topic :: Package Development :: Build Tools',
	'License :: OSI Approved :: MIT License',
	'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    ],
	packages = find_packages(),
	install_requires = ["Orange-Bioinformatics"],
	include_package_data = True,
)

