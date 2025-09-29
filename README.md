# Variable-DNA-regions-Identification

Program created to identify variable DNA region in a species/strain and especially indels (insertion and deletion) from analysis of different genomes. Mainly used to implement a new MLVA strategy.
<br>This program has been created using *Cryptosporidium hominis* genomes.

## How to use it

To use this program you have to create a folder and put all the genomes you want to be analysed inside. The genomes files must be pileup or mpileup type.

## Functions overview

The first and second step are there to determined the variable regions of the genome. Their output is the pileup and csv files 'analysis'.
<br>The third step randomizes and chooses a number of windows then classify the samples using only the informations contained in the choosen windows. This can be use to implement a MLVA strategy.

#### `main`
Calls all the functions and manages the folder and files.

#### `log`
Display the time for each step so following the program become easier. 


### **First step** : to find and regroup the indels

#### `researche_change`
Will find and count all the indels, then finds the most recurring one. Also cleans the raw pileup file.

#### `regroup_mut`
Will regroup all the indels found for one window together, for every windows.

#### `process_sample`
Calls researche_change & regroup_mut for every file

#### `find_window`
Creates a dictionnary for every windows in every sample & writes all the windows that contain indels


### **Second step** : to calculate the score and to create the output

#### `calc_score`
Create a score file in pileup for each window, with their position, their score and all the indels found for every sample

#### `visual`
Creates the first csv file "analysis" from the pileup generated earlier (seeing the results as a human become easier)

 
### **Third step** : to randomize

#### `group_numpy`
Classify the samples with indels found

#### `random`
Select random windows and send them to group_numpy to classify them

#### `main`

some parts


## Author

Created by Lou PLANTEROSE and Romain COPPÉE.  
GitHub: [Lou PLANTEROSE](https://github.com/Looplye)
