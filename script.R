
# General Information:
# ====================
#
# For a more comprehensive introduction to this script, see 'Readme.md' on:
# https://github.com/BalduinLandolt/createStemmaR
#
# NB:
# ---
# Every now and again, phangorn::bootstrap.phyDat() seems to cause crashes.
# So far, I found no way of reproducing that bug reliably.
# Since it sporadicly happens with the pre-prepared data, which does not change,
# it can hardly be my fault. (Persumably it is caused by the random values
# that are used by the bootstrap algorithm.)
# Should this happen, please just restart RStudio, set the working directory again,
# and re-run the script. Then it should work.


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
# (method that does that is specified in functions.R)
trees = make_trees(data, name = "trivial")
# and plot the according dendrograms
# (method that does that is specified in functions.R)
plot_trees(trees)





# generate a random sample of data
# (you'll get complete nonsense, of course.)
# (method that does that is specified in functions.R)
data = get_random_data()

# print head of data
head(data)

# calculate and plot the same dendrograms from this data
trees = make_trees(data, name = "random")
plot_trees(trees)




