# Variable-DNA-regions-Identification

Program created to identify variable region in a genome and especially indels (insertion and deletion). Mainly used to implement a new MLVA strategy.
This program has been created using *Cryptosporidium hominis* genomes.

## How to use it

To use this program you have to create a folder and put all the genomes you want to be analysed inside. The genomes files must be pileup or mpileup type.

## Functions overview

### First step : to find and regroup the indels

#### `researche_change`
Will find and count all the indels, then finds the most recurring one. Also cleans the raw pileup file.

#### `regroup_mut`
Will regroup all the indels found for one window together, for every windows.

#### `process_sample`
Calls researche_change & regroup_mut for every file

#### `find_window`
Creates a dictionnary for every windows in every sample & writes all the windows that contain indels


### Second step : to calculate the score and to create the output

#### `calc_score`
Create score file

#### `visual`
Creates first csv file "analysis"

 

#Classify the samples with indels found
def group_numpy(call_matrix):


#Select random windows and send them to group_numpy to classify them
def random

## Author

Created by Lou PLANTEROSE and Romain COPPÉE.  
GitHub: [Lou PLANTEROSE](https://github.com/Looplye)
