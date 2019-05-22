
# General Information:
# ====================
#
# For a more comprehensive introduction to this script, see 'Readme.md' on:
# https://github.com/BalduinLandolt/createStemmaR
#


# Run instructions:
# =================
#
# - set working directory to where the data lies, or adjust paths accordingly.
# - install required packages
# - run script step by step


# imports:
# ========
#
# Required packages:
# - ape
#   [run: install.packages("ape") ]
# - phangorn
#   [run: install.packages("phangorn") ]





# Script starts here
# ==================




# load packages
library(ape)
library(phangorn)

# load functions
source("functions.R")




# load pre-made trivial data sample
data = read.csv("data/trivial.csv", header = T, sep = ";", stringsAsFactors = F, skip = 0)

# print head of data
head(data)

# from that data, calculate trees using several different algorythms,
trees = make_trees(data, name = "trivial")
# and plot the according dendrograms
# (method that does that is specified below)
plot_trees(trees)





# generate a random sample of data
# (you'll get complete nonsense, of course.)
data = get_random_data()

# print head of data
head(data)

# calculate and plot the same dendrograms from this data
trees = make_trees(data, name = "trivial")
plot_trees(trees)




