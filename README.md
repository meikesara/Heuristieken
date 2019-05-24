# Protein Pow(d)er by 'EiwitVouwers'
#### Definition of protein -- TODO


Hier staat een korte beschrijving van het probleem evt. met plaatje.
(linken naar de site)

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
* Will only run in respectable time for proteins with <b>20 amino acids or less</b>


##### Random
To run the random algorithm use (protein is a string of amino acids):
```
python main.py random protein
```
Thereafter, you are asked for the dimension in which the folding should be performed (either 2D or 3D) and the amount of runnings that should be performed.

If one running is performed, the folded protein will be visualized and an estimate of the best stability will also be given.</br>
If multiple runnings are performed, you will be prompted that the results will be saved in files. In ```proteinRandom.txt``` the coordinates on which the amino acids are placed for each protein found are saved. In ```stabilityRandom.txt``` the stability for each running is saved.</br> Thereafter, both the best found stability and an estimate of the best stability are returned.


##### Hill Climber
To run the hill climber algorithm use (protein is a string of amino acids):
```
python main.py hillclimber protein
```
Thereafter, you are asked for the dimension in which the folding should be performed (either 2D or 3D), the amount of runnings that should be performed, and the amount of iterations per running that should be performed.

If one running is performed, the folded protein and the change of the stability over the iterations will be visualized and the estimate of the stability will also be given.</br>
If multiple runnings are performed, you will be prompted that the results will be saved in files. In ```proteinHillClimber.txt``` the coordinates on which the amino acids are placed (for the best folded protein of each run) are saved. In ```stabilityHillClimber.txt``` the best found stability for each running is saved.</br>
Thereafter, both the overall best found stability and an estimate of the best
possible stability are returned.

##### Simulated Annealing
To run the simulated annealing algorithm use (protein is a string of amino acids):
```
python main.py simulated protein
```
This is a logarithmic simulated annealing algorithm and the cooling schedule used is: ```D/ln(iterations + 2) - D/ln(10^6)```

You are asked for the dimension in which the folding should be performed (either 2D or 3D), the amount of runnings that should be performed, the amount of iterations per running that should be performed, and the D (as in the function for the cooling schedule above).

If one running is performed, the folded protein and the change of the stability, temperature, and acceptance rate over the iterations will be visualized and the estimate of the stability will also be given.</br>
If multiple runnings are performed, you will be prompted that the results will be saved in files. In ```proteinHillClimber.txt``` the coordinates on which the amino acids are placed (for the best folded protein of each run) are saved. In ```stabilityHillClimber.txt``` the best found stability for each running is saved.</br>
Thereafter, both the overall best found stability and the estimated stability are returned.


## Authors

* Nicole Jansen
* Meike Kortleve


## Acknowledgments

* StackOverflow
* minor programmeren van de UvA
* TODO -- referentie Lesh
