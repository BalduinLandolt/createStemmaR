# createStemmaR

Testing possibilities of generating phylogenetic stemmata with R.


### Run instructions

- Set Working directory to where the script path. All input data needs to be placed in a folder `data` in the same directory. (Otherwise, paths in the script need to be adjusted accordingly.)
- Install required packages.


### Required Packages

- ape
- phangorn


### What it does

The script reads a trivial sample data set. (i.e. aligned data matrix)

From the input data it calculates an NJ-tree and plots it.


### To Do

- Random data
- more kinds of tree calculation


## Workflow

The following workflow is intended for this script.  
(For phylogenetic analysis of humanities data. E.g. Folk Tales or Manuscripts.)

1. Align data in Excel Sheet:  
Column 1: underscore plus unique Name of data row (manuscript, text, ...); name should be exactly 10 characters long.  
Row 1: Titles for columns (for complex data, like texts, this could be describing, what the differences refer to.)

2. Analyze data:  
Under every row of data, add an extra empty row.  
Add the same unique name, just without the underscore.  
Then, for every column, analyze the data. Assign values for data; use digits 0-9 for values, "-" for missing element; "?" for unclear element.

3. Extract data:  
Sort table alphabetically. Thanks to the underscore, the taxonomy gets separated from the original data.   
Export the taxonomy as .csv file.

4. Run script.R  
The script should reate a valid .nex file of this alignment, calculate the matching tree, save the tree, and plot it.


