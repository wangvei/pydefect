# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Automatic constructor of point-defect enveironment, and analyzer of vasp 
  calculation results.
* 0.2.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

#### Installation using pip

After cloning the repository, it is possible to install `pydefect` using the python package manager pip.
To do so, 
run this command in the directory containing setup.py:

`pip install ./`

This will install `pydefect` and any dependant packages
not already have installed on your system. 

Note that the `scikit-image` package is required to locate interstitials but can be omitted if this feature is not required.

Sometimes errors will be given if a specific version of a package is not 
installed. In this case try installing the exact versions that have been 
tested for use with `pydefect` with the following command:

`pip install -r requirements.txt ./`

This will install the package versions listed in the requirements.txt file.
To prevent interference with other programmes, it is advised that a package 
management system like 
[conda](https://docs.conda.io/projects/conda/en/latest/index.html) is used. 
This allows you to install dependancies in a particular environment so that 
they can be managed and recorded more easily. The commands above can then be 
executed in a conda environment (after installing pip in that environment).

For more information on how dependancies are managed in this branch see this [blog post](https://medium.com/@boscacci/why-and-how-to-make-a-requirements-txt-f329c685181e).



* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
