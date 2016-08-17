#!/usr/bin/env python

from setuptools import setup, find_packages
from distutils.core import setup


setup(
	name = "ymap",
	
	version = "2.0.1",
	
	description = "An automated method to map yeast variants to protein modifications and functional regions",
	
	url = "https://github.com/CSB-KUL/yMap",
	
	author = "Ahmed Arslan and Vera van Noort, KU Leuven, Belgium",
	
	author_email = "vera.vannoort@kuleuven.be",
	
	maintainer = "Ahmed Arslan",
	
	maintainer_email = "ahmed.arslan@kuleuven.be",
	
	license = "MIT",
	
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

    	Package_data = {},

	packages = ['ymap'],

	install_requires = ["Orange-Bioinformatics", "Orange"],
	
	include_package_data = True,
	
	entry_points={
        	'console_scripts': [
            
            		'ydata=ymap.ymap:data',
            
            		'ygenes=ymap.ymap:ymap_genes',
            
        		 'yproteins=ymap.ymap:ymap_proteins',
            
            		'yweb=ymap:web'

        ],
    }
)


