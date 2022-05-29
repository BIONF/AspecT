# XspecT - Acinetobacter Species Assignment Tool
<p align="center">
<img src="https://github.com/BIONF/XspecT/blob/main/static/Logo.png" height="50%" width="50%">
</p>
XspecT is a Python-based tool to taxonomically classify [Acinetobacter](https://en.wikipedia.org/wiki/Acinetobacter) sequence-reads (or assembled genomes) on the species and/or sub-type level using [Bloom Filters](https://en.wikipedia.org/wiki/Bloom_filter) and a [Support Vector Machine](https://en.wikipedia.org/wiki/Support-vector_machine). It also identifies existing [blaOxa-genes](https://en.wikipedia.org/wiki/Beta-lactamase#OXA_beta-lactamases_(class_D)) and provides a list of relevant research paper for further information.
<br/><br/>

XspecT utilizes the uniqueness of kmers and compares extracted kmers from the input-data to a reference database. Bloom Filter ensure a fast lookup in this process. For a final prediction the results are classified using a Support Vector Machine. 
<br/>

Local extensions of the reference database is supported.
<br/>

The tool is available as a web-based application and a smaller command line interface.

## Table of Content
- [Installation and Usage](#installation-and-usage)
    + [Python Modules - Install Requirements](#python-modules---install-requirements)
    + [List of used Modules for Python (3.10):](#list-of-used-modules-for-python--310--)
    + [How to run the web-app: Local Deployment](#how-to-run-the-web-app--local-deployment)
- [Input Data](#input-data)
- [Walkthrough](#walkthrough)
- [Contributors](#contributors)
- [About this project](#about-this-project)


## Installation and Usage
#### Python Modules - Install Requirements
On Linux you need the python-dev package:
```
sudo apt install python3.10-dev
```
XspecT requires the latest 64 bit Python version and a list of Python Modules (see below).
```
pip install -r requirements.txt
```
#### List of used Modules for Python (3.10):
- Flask
- Flask-Bcrypt
- Flask-Login
- Flask-WTF
- WTForms
- Werkzeug
- Bcrypt
- Biopython
- bitarray
- mmh3
- numpy
- pandas
- requests
- scikit-learn

#### How to run the web-app: Local Deployment
Run the following command lines in a console, a browser window will open automatically after the application is fully loaded.

MAC/Linux:
```
$ export FLASK_APP=flaskr
$ export FLASK_ENV=development
$ python app.py
```
Windows cmd:
```
set FLASK_APP=flaskr
set FLASK_ENV=development
python app.py
```
To start the command line interface of XspecT use the following command line (on all OS):
```
Python XspecT_mini.py
```

## Input Data
XspecT is able to use either raw sequence-reads (FASTQ-format .fq/.fastq) or already assembled genomes (FASTA-format .fasta/.fna). Using sequence-reads saves up the assembly process but high-quality reads with a low error-rate are needed (e.g. Illumina-reads).

The amount of reads that will be used has to be set by the user when using sequence-reads. The minimum amount is 5000 reads for species classification and 500 reads for sub-type classification. The maximum number of reads is limited by the browser and is usually around ~8 million reads. Using more reads will lead to a increased runtime (xsec./1mio reads).

## Walkthrough
A detailed walkthrough with examples is provided in Xspect’s wiki.

## Contributors
- Sam Gimbel
- Bardja Djahinschiri
- Vinh Tran
- Ingo Ebersberger

## About this project
This project is an attempt to support hospital staff in a possible A. baumannii outbreak. A. baumannii can build up antibiotic resistance and can cause deadly nosocomial infections.
This is a bachelor thesis project; no warranty is given. Check the license for more information.
constructive criticism/feedback always welcomed!
  
