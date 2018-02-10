from setuptools import setup, find_packages

setup(
    name='PyDefect',
    version='0.1.dev0',
    author='Yu Kumagai',
    author_email='yuuukuma@gmail.com',
    url='https://scholar.google.co.jp/citations?user=xST4MSEAAAAJ&hl=ja',
    packages=['bin','test',],
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ], install_requires=['numpy', 'pymatgen']
)
