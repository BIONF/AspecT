# ClAssT - Clone-type Assignment Tool for A.baumannii
# AspecT - Acinetobacter Species Assignment Tool
ClAssT and AspecT are both using Bloom Filters and Support Vector Machines to classify sequence-reads (.fq-files) or assembled genomes (.fasta or .fna files).
ClAssT assigns the Input-Data to one of the eight international clones (IC) of A.baumannii, if the input sequence is part of any of the eight clones.
AspecT assigns the Input-Data to one of the 82 Acinetobacter Species, if the Input-Data is similar enough to one of the Reference-BloomFilter.
The Tool is web-based.

# How ClAssT works

![alt text](https://github.com/Dominik0304/AspecT/blob/master/static/Workflow_ClAssT.png)

The Tool uses Bloom Filters to store IC-specific reference k-meres. The k-meres of the input-sequences will be checked for membership in all of the selected IC's. Hits are counted and then divided by the number of total tested k-meres of the input-sequence. This produces a 'Score-Vector' with values between 0 and 1. This Vector will then be classified by Support Vector Machines (SVM).

If you are using sequence-reads for this analysis, then you can avoid the assembly-process. On the other hand, this Tool need high quality sequence-reads because its using exact matches of k-meres.

The requirements are down below and in the requirements.txt file.

For security reasons, you need to set up an individual security key and a username/password for a login. More information in the following instructions.

# How to use ClAssT
## Assigning Files
  <b> 1) Choose a file and submit</b><br>
  <b> 2) wait </b><br>
  <b> 3) get Results </b><br>
<p align="center">
  <img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/How2Use.png" height="50%" width="50%">
</p>

## Modify the Tool
Searchable Filters can be added or removed in the 'Export Options'- Section on the website. Login to that area (see section 'Change Password and Username' in this readme) and follow the upcoming steps.

### Adding Genomes
1) Choose a name
2) Select the file that contains the Genome from your Computer
3) Add/Expand the SVM Training-Vectors: <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 3.1) Change the 'Score_new' to the corresponding Value between 0 and 1 <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 3.2) Add the new Score-Vectors of the genome with the corresponding Value between 0 and 1. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The format must be ['Filename', Score-IC1, Score-IC2,..., Score_of_new, Label] like all other Vectors showen below.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The Label in the last column must match the previous entered name in step 1). <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; There needs to be at least one Training-Vector for the new Genome. <br>
<p align="center">
  <img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/AddFilter.png" height="50%" width="50%">
</p>

### Removing Genomes
All deletable Genomes are shown. Copy the name of the one you want to delete, paste in into the text field and submit.

### Add OXA-Genes
Adding OXA-Genes is almost similar to adding Genomes, but without any Score-Vectors. Add the .fasta-file with the gene, put a name into the text-field and submit.

### Remove OXA-Genes
All deletable OXA-Genes are shown. Copy the name of the one you want to delete, paste in into the text field and submit.

### Modify Trainingdata for the SVM
This function allows you to change the Trainingvectors for the SVM. There must be at least one vector per label!
<p align="center">
  <img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/modify_vecs.png" height="50%" width="50%">
</p>

# How AspecT works

AspecT works roughly the same like ClAssT. For a Taxonomic Assignment to be performed on the species level the Input-Data is processed with the same File-Reader and the Bloom Filter are searched through the same modules. For each species up to 4 different assembled genomes (if availabe) were used as reference-data for the Bloom Filters.
Only a small amount of k-mers are needed for the species assignment. AspecT uses only every 500th k-mere of a assembled genome and only every 10th k-mer of a Sequence-Read.
The Support Vector Machine (SVM) was further adjusted using the radial basis function as kernel-function with a regularization parameter of C = 1.5

# How to use AspecT
## Assigning Files

<b> 1) Choose a .fna/.fasta/.fastq-file and submit</b><br>
<b> 2) wait (this should take just a few seconds)</b><br>
<b> 3) get Results </b><br>
<p align="center">
<img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/HowtouseAspecT.png" height="50%" width="50%">
</p>

# Add new Acinetobacter Species
If new Species are discovered it is possible to add them to AspecT.

## Adding species
AspecT needs Reference-Data for the BloomFilter. For an accurate Assignment 4 Assemblys for each new species are recommended.
Those Assemblys need to be concatenated into one File and should not contain more than ~11.9 Million distinct kmers. You can use the Tool [Jellyfish](https://github.com/gmarcais/Jellyfish) to count kmers.

1) Make sure that the file-name is the species-name
2) Copy and Paste the (concatenated-)file in the folder filter/new_species
<p align="center">
<img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/AddSpecies1.png" height="50%" width="50%">
</p>



## Adding training-data
AspecT needs Training-data (Genome Assemblies) for each supported species. If no separate Training-data is available you can use the BloomFilter-Reference-Data.

1) Make sure that the file-name contains an Assembly-Accession (for better read-ability in the csv-file)
2) Copy and Paste the file in the folder Training_data/genomes
<p align="center">
<img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/AddSpecies2.png" height="50%" width="50%">
</p>

## Run the Script

1) Run the Script Add_Species.py with the Species-Names as additional Parameters
```
python Add_Species NewSpecies1 NewSpecies2
```
2) The new added Species will be trained into a BloomFilter
3) New SVM-Training-Data will be generated
4) After the Script is finished make sure to delete the Assembly-Files from the filter/new_species Folder
5) Restart the Tool to apply the new changes


# Setup
### Python Modules - Install Requirements
Using a virtual environment is recommended.
```
pip install -r requirements.txt
```

#### List of used Modules for Python (3.8)
Flask	1.1.2

Flask-Bcrypt	0.7.1

Flask-Login	0.5.0

Flask-WTF	0.14.3

WTForms	2.3.1

Werkzeug	1.0.1

bcrypt	3.1.7

biopython	1.76

bitarray	1.2.1

mmh3	2.5.1

numpy	1.18.2

pandas	1.0.3

requests	2.23.0

scikit-learn	0.23.1

### Setup

#### Change Secret Key
Because of security reasons, you need to give this Tool a new Secret Key. Change the 'change_me' in the settings.cfg in the config folder to your own Secret Key :
<p align="center">
  <img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/secretkey.png" height="30%" width="30%">
</p>

You can generate a random Secret Key of variable length by using Python:
```
>>> import os
>>> os.urandom(20)
b'\x80\xf8\xfe\xbe\xb5t*{\x88\xdc\xb3z\x17\xacz\xeasM\xf7\xd4'
```
#### Change Password and Username in 'Expert Options'
To gain or to limit the access to the 'Expert Options' you need to change the Username and Password by using the change_password.py script:
<p align="center">
  <img src="https://github.com/Dominik0304/AspecT_ClAssT/blob/main/Instructions/pictures/change_pw.png" height="50%" width="50%">
</p>

# How to run the App: Local Deployment

## MAC/Linux:

```
$ export FLASK_APP=flaskr

$ export FLASK_ENV=development

$ python app.py
```

## Windows cmd:

```
set FLASK_APP=flaskr

set FLASK_ENV=development

python app.py
```

## More about the Deployment

[Documentation](https://flask.palletsprojects.com/en/master/tutorial/factory/)

# Server Deployment

[Server Deployment](https://flask.palletsprojects.com/en/master/deploying/)

[Host a app for free by using Heroku](https://www.heroku.com/): NOT RECOMMENDED (yet) - Heroku is an easy (and free) method to host Flask Apps, but Heroku will throw [Request timeout errors](https://devcenter.heroku.com/articles/error-codes#h12-request-timeout) if the server takes longer than 30 seconds to respond. Some Options take longer than 30 seconds, therefore hosting on Heroku is problematic in this case. A job scheduler can solve this problem and is planned for the future.

# Web-Interface
The code for the Web-interface was taken from [Corey M. Schafer](https://github.com/CoreyMSchafer/code_snippets/tree/master/Python/Flask_Blog). Some modifications have been made, sources can be found in the code.  

# About this project
This project is an attempt to support hospital staff in a possible A.baumannii outbreak. A.baumannii is able to build up antibiotic resistance and can cause deadly nosocomial infections.

This is a bachelor thesis project, no warranty is given. Check the licence for more information.

constructive criticism/feedback always welcomed!
"# AspecT_ClAssT" 
