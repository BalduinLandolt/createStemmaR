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

For both of them, it calculates a couple of different trees, using multiple alorithms.

Then it plots dendrograms for these trees.


### To Do

So far so well.


## Workflow

The following workflow is intended for this script.  
(For phylogenetic analysis of humanities data. E.g. Folk Tales or Manuscripts.)

1. Align data in Excel Sheet:  
Column 1: underscore plus unique Name of data row (manuscript, text, ...); name should be exactly 10 characters long. The following columns contain the actual data.  
Row 1: Titles for columns (for complex data, like texts, this could be describing, what the differences refer to). Rows 2 onward contain data from one artifact (text/manuscript).

2. Analyze data:  
Under every row of data, add an extra empty row.  
Add the same unique name, just without the underscore.  
Then, for every column, analyze the data. Assign values for data; use digits 0-9 for values, "-" for missing element; "?" for unclear element.

3. Extract data:  
Sort table alphabetically. Thanks to the underscore, the taxonomy gets separated from the original data.   
Export the taxonomy as .csv file.

4. Run script.R  
The script should create a valid .nex file of this alignment, calculate the matching tree, save the tree, and plot it.  
For this to happen, you'd need the same method calls as are used for the trivial or the random data (i.e. `read.csv()`, `make_trees()`, `plot_trees()`).
