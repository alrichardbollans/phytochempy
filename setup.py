from setuptools import setup, find_packages

setup(
    name='phytochempy',
    version='1.1',
    packages=find_packages(),
    package_data={"phytochempy": ["compound_properties/inputs/*"]},
    install_requires=[
        'pandas',
        'numpy',
        'tqdm',
        'rdkit',
        'standardiser',
        'wcvpy>=1.3.2'
    ],
    extras_require={'knapsack': ["html5lib", 'beautifulsoup4', 'cirpy'],
                    'compound_metrics': ['rdkit', 'chembl_webresource_client'],
                    },
    url='https://github.com/alrichardbollans/phytochempy',
    license='GNU v.3',
    author='Adam Richard-Bollans',
    description='A python package to download and analyse phytochemical data ',
    long_description=open('README.md', encoding="utf8").read()
)
