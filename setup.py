#!/usr/bin/env python

import sys
import os
import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
from distutils.core import setup

if sys.version_info >= (3,1):
	install_require = ["Orange3"]

elif sys.version_info < (3,0):
	install_require = ["Orange"]


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
	name = "ymap",
	
	version = "2.0.6a",
	
	description = "An automated method to map yeast variants to protein modifications and functional regions",
	
	url = "https://github.com/CSB-KUL/yMap",
	
	author = "Ahmed Arslan and Vera van Noort, KU Leuven, Belgium",
	
	author_email = "vera.vannoort@kuleuven.be",
	
	maintainer = "Ahmed Arslan",
	
	maintainer_email = "ahmed.arslan@kuleuven.be",
	
	license = "MIT",
	
	long_description=read('README.md'),
	
	keywords = (
		
		'Biofinformatics',
		
		'Proteins variants',
		
		'PostTransational-Modifications (PTMs)',
		
		'Proteins domains',
		
		'Secondary Structures',
		
		'Protein-Protein Interactions',
		
		'Genomics',
		
		'Proteomics',
		
		'Gene Ontology',
		
		'yMap',
	
		),

	classifiers=[ 
	'Intended Audience :: Science/Research',
	'Topic :: Education',
	'License :: OSI Approved :: MIT License',
	'Programming Language :: Python :: 2.7',
    	'Programming Language :: Python :: 3',
    	'Programming Language :: Python :: 3.3',
    	'Programming Language :: Python :: 3.4',
    	'Programming Language :: Python :: 3.5',
    ],

    	install_requires = [install_require, "Orange-Bioinformatics"],

	packages = ['ymap'],
	
	include_package_data = True,
	
	zip_safe = False,
	
	entry_points={
        	'console_scripts': [
            
            		'ydata = ymap.ymap:data',
            
            		'ygenes = ymap.ymap:ymap_genes',
            
        		 'yproteins = ymap.ymap:ymap_proteins',
            
            		'yweb = ymap.ymap:web',
            		
            		'ytest = ymap.tests:main',

        ],
    }
)


