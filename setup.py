from setuptools import setup

###
# @author: Pietro Cinaglia
# @mail: cinaglia (at) unicz (dot) it
# @github: https://github.com/pietrocinaglia/corte
###

setup(
    name='CoRTE',
    version='0.1',    
    description='An open source and user-friendly tool for COnstructing Real-world TEmporal networks from genotype-tissue expression data (CoRTE)',
    url='https://github.com/pietrocinaglia/corte',
    author='Pietro Cinaglia',
    author_email='cinaglia@unicz.it',
    license='MIT',
    packages=['corte'],
    install_requires=['matplotlib>=3.7.2',
                      'networkx>=2.8.6',
                      'pandas>=2.0.3',
                      'Requests>=2.32.3'
                      'scipy>=1.14.1',          
                      ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.8',
    ],
)