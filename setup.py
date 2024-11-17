from setuptools import setup, find_packages

setup(
    name='GISTA',
    version='1.0',
    packages=find_packages(),  # Automatically find the 'GISTA' package
    install_requires=[
        'Click', 'scipy', 'numpy', 'matplotlib', 'pandas', 'seaborn', 'tqdm', 'logomaker'
    ],
    entry_points='''
        [console_scripts]
        GISTA=GISTA:cli
    '''
)