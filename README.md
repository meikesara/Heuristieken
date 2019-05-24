# Protein Pow(d)er by 'EiwitVouwers'
#### Definition of protein
TODO
Hier staat een korte beschrijving van het probleem evt. met plaatje.

## Getting Started

### Prerequisites

This codebase has been written in [Python3.7.1 and Python3.7.2](https://www.python.org/downloads/). The file requirements.txt contains all the necessary packages to successfully run this code. These can be installed using pip with the following instruction:

```
pip install -r requirements.txt
```

### Structure

Alle Python scripts staan in de folder Code. In de map Data zitten alle input waardes en in de map resultaten worden alle resultaten opgeslagen door de code.

### Usage
##### Constructive
To run the constructive (depth first) algorithm use (protein is a string of amino acids):

```
python main.py constructive protein
```
Example:
```
python main.py constructive HHPHHHPH
```
Returns the optimal stability.

<i>Notes:</i>
* This algorithm only works in 2D
* Will only run in respectable time for proteins with <b>20 amino acids or less<b/>


##### Random
To run the random algorithm use (protein is a string of amino acids):
```
python main.py random protein
```
Thereafter, you are asked for the dimension in which the folding should be performed (either 2D or 3D) and the amount of runnings that should be performed. If one 

Returns


##### Hill Climber
To run the hill climber algorithm use (protein is a string of amino acids):
```
python main.py hillclimber protein
```

## Authors

* Nicole Jansen
* Meike Kortleve


## Acknowledgments

* StackOverflow
* minor programmeren van de UvA
