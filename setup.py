from setuptools import setup, find_packages

setup(
    name='phytochempy',
    version='0.1',
    packages=find_packages(include=['phytochempy'], exclude=['library_info_and_data_import']),
    package_data={"phytochempy": ["compound_properties/inputs/*"]},
    install_requires=[
        'pandas>=2.1.4',
        'numpy',
        'urllib',
        'cirpy'
    ],
    extras_require={'knapsack': ["html5lib", 'beautifulsoup4','tqdm'],
                    'metabolite_metrics': ['rdkit', 'chembl_webresource_client'],
                    },
    url='https://github.com/alrichardbollans/PhytoChemicalDiversity',
    license='GNU v.3',
    author='Adam Richard-Bollans',
    description='A package for analysing chemical diversity in vascular plants',
    long_description=open('README.md', encoding="utf8").read()
)
