# createStemmaR
(By Balduin Landolt)

This project aims to test possibilities of generating phylogenetic stemmata with R.

### Project Structure
The expected project structure must be:

`root/script.R`  
`root/functions.R`  
`root/data/trivial.csv`

### Run instructions

- Ensure that the project structure is correct and all files are where they are expected to be.  
Easiest is to clone the git repository.  
If another structure is wanted, all file paths in the scripts need to be adjusted accordingly.

- Install required packages.

- Run `script.R`



### Required Packages

- ape
- phangorn


### What it does

The script reads a provided trivial sample data set and creates a random one.

For both of them, it calculates a couple of different trees, using multiple algorithms.

Then it plots dendrograms for these trees.


### To Do

So far so well.


## Workflow

The following workflow is intended for this script.  
(For phylogenetic analysis of humanities data. E.g. Folk Tales or Manuscripts.)

1. Align data in Excel Sheet:  
Column 1: underscore plus unique Name of data row (manuscript, text, ...); name should be exactly 10 characters long. The following columns contain the actual data.  
Row 1: Titles for columns (for complex data, like texts, this could be describing, what the differences refer to). Rows 2 onward contain data from one artefact (text/manuscript).

2. Analyse data:  
Under every row of data, add an extra empty row.  
Add the same unique name, just without the underscore.  
Then, for every column, analyse the data. Assign values for data; use digits 0-9 for values, "-" for missing element; "?" for unclear element.

3. Extract data:  
Sort table alphabetically. Thanks to the underscore, the taxonomy gets separated from the original data.   
Export the taxonomy as .csv file.

4. Run script.R  
The script should create a valid .nex file of this alignment, calculate the matching tree, save the tree, and plot it.  
For this to happen, you'd need the same method calls as are used for the trivial or the random data (i.e. `read.csv()`, `make_trees()`, `plot_trees()`).



## Background

The idea of using biological metaphors to describe patterns in the humanities is by no means new.  
Since the 19th century, languages have understood as organic, living beings; their relationships were described like genealogic relations. Similarly, in textual criticism, the relationship between manuscripts were seen as genealogical - most famously by the German scholar Karl Lachmann; although he never published a _stemma_ himself (the earliest _stemma_ published, to my knowledge, is the one by the Swedish scholars C.J. Schlyter and H.S. Collin in their edition of _Westgöta-Lagen_ from 1827), the _stemma codicum_ is based on the Lachmannian concept of genealogical manuscript relationship. It applies the visualisation of the family tree to manuscripts.  
Lotman's application of the term 'morphology' to fairy tails is a borrowing from biology as well. And with the rise of _oral theory_ in the middle of the middle of the 20th century, the categorization of folk literature was conducted under the genealogical premises.

It is no surprise that techniques from bioinformatics have been adopted for cultural analysis; especially the method of computationally creating _stemmata_ from aligned taxa, has proved very promising both in manuscript studies and in folklore studies.  
The fact that these methods usually create so called 'unrooted' dendrograms, is very much in line with modern approaches like material philology, that have previously rejected the concept of the Lachmannian _stemma_ for being too hierarchical, and thus - more or less subconsciously - reinforcing the rejected concept of the archetype as the ultimate goal of textual criticism.  
(Or, comparably, for example Dumézil's Indo-European mythographical theories in the studies of orally transmitted texts.)

Seeing the promising results of these approaches, it is a pity to see how few such studies have actually been conducted. That is, for once, due to the fact how much work needs to be invested, to convert a human made artefact into a machine readable data alignment.  
Another may well be that it's not exactly easy for scholars of humanities to apply these methods in a very practical sense: There are numerous tools that can do this (e.g. the Phylip package, PAUP or SplitsTrees) but they are not all really user friendly (e.g. Phylip being a package of numerous command-line based tools, first released in 1980). The tools usually support an overwhelming number of algorithms with even more variable parameters.  
For a first time user, these tools are complete black boxes, because the documentations are rather meager, and definitely not targeted at users who are not trained professionals in the field of bioinformatics.  
Also, the general problem with using prefabricated tools applies here: By using such a tool, one has to rely on it being available, on it working correctly, being stable, and so forth. And using a tool means giving up control over the data - thus, even after consulting all available documentation, a certain extend of black box remains.

For all of these reasons, and because I intend to work with phylogenetic methods on Scandinavian and European folk ballads in the near future, I decided to test out the possibilities of calculating and plotting dendrograms using the programming language R.
