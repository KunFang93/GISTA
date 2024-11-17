from setuptools import setup

setup(
    name='GISTA',
    version='1.0',
    py_modules=['GISTA'],
    install_requires=[
        'Click','scipy','numpy','matplotlib','pandas','seaborn','tqdm','logomaker'
    ],
    entry_points='''
    [console_scripts]
    GISTA=GISTA:cli
    '''
)