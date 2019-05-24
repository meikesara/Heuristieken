# Protein Pow(d)er by 'EiwitVouwers'
Proteins are long strings of amino acids. The correct folding of a protein is crucial for its functioning. If a protein is not folded correctly, it could lead to pathologies. It is therefore important to know in which way the protein could be correctly folded, especially for fields such as pharmaceutics.

The problem we had to solve was therefore to find a solution with which proteins could be folded as efficiently as possible. The whole problem is described [here][1] (in Dutch). A quick summary: we have to find a solution for the Hydrophobic-Polar protein folding model. In this model, all amino acids are classified as either hydrophobic (H) or polar (P). Additionally, when two hydrophobic amino acids lay adjacent, they will form a bond, which decreases the stability by 1. Here, if the stability is lower, it means the protein is folded more stably. Additionally, we have to add cysteine (C) amino acids, and a bond between two cysteine amino acids will decrease the stability by 5. We only have to fold the proteins on a square 2D grid and a cubic 3D grid and with angles of only 90 degrees.

[1]: https://heuristieken.nl/wiki/index.php?title=Protein_Pow(d)er

## Getting Started

### Prerequisites
This codebase has been written in [Python3.7.1 and Python3.7.2](https://www.python.org/downloads/) in Windows. The file ```requirements.txt``` contains all the necessary packages to successfully run this code. These can be installed using pip with the following instruction:

```
pip install -r requirements.txt
```

### Structure
In the folder "algorithms" are the Python scripts with the different algorithms that could be used. In the folder "classes" are the Python scripts with the classes for amino acids and proteins. In the folder "helper" are the Python scripts with helper functions such as a visualizer and code for estimating the lower bound of the stability.

### Usage
#### Constructive
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


#### Random
To run the random algorithm use (protein is a string of amino acids):
```
python main.py random protein
```
Thereafter, you are asked for the dimension in which the folding should be performed (either 2D or 3D) and the amount of runnings that should be performed.

If one running is performed, the folded protein will be visualized and an estimate of the best stability will also be given.</br>
If multiple runnings are performed, you will be prompted that the results will be saved in files. In ```proteinRandom.txt``` the coordinates on which the amino acids are placed for each protein found are saved. In ```stabilityRandom.txt``` the stability for each running is saved.</br> Thereafter, both the best found stability and an estimate of the lower bound of the stability are returned.


#### Hill Climber
To run the hill climber algorithm use (protein is a string of amino acids):
```
python main.py hillclimber protein
```
Thereafter, you are asked for the dimension in which the folding should be performed (either 2D or 3D), the amount of runnings that should be performed, and the amount of iterations per running that should be performed.

If one running is performed, the folded protein and the change of the stability over the iterations will be visualized and the estimate of the lowe bound of the stability will also be given.</br>
If multiple runnings are performed, you will be prompted that the results will be saved in files, which could be found in the folder results. In ```proteinHillClimber.txt``` the coordinates on which the amino acids are placed (for the best folded protein of each run) are saved. In ```stabilityHillClimber.txt``` the best found stability for each running is saved. Thereafter, both the overall best found stability and an estimate of the lower bound of the stability are returned.

#### Simulated Annealing
To run the simulated annealing algorithm use (protein is a string of amino acids):
```
python main.py simulated protein
```
This is a logarithmic simulated annealing algorithm and the cooling schedule used is: ```D/ln(iteration + 2) - D/ln(10^6)```. If an amount of iterations greater than 10^6 - 2 is chosen, the temperature will be: ```D/ln(10^6 - 1) - D/ln(10^6)```.

You are asked for the dimension in which the folding should be performed (either 2D or 3D), the amount of runnings that should be performed, the amount of iterations per running that should be performed, and the D (as in the function for the cooling schedule above).

If one running is performed, the folded protein and the change of the stability, temperature, and acceptance rate over the iterations will be visualized and the estimate of the stability will also be given.</br>
If multiple runnings are performed, you will be prompted that the results will be saved in files, which could be found in the folder results. In ```proteinSimulatedAnnealing.txt``` the coordinates on which the amino acids are placed (for the best folded protein of each run) are saved. In ```stabilitySimulatedAnnealing.txt``` the best found stability for each running is saved. Thereafter, both the overall best found stability and an estimate of the lower bound of the stability are returned.


## Authors

* Nicole Jansen
* Meike Kortleve


## Acknowledgments

* StackOverflow
* Minor Programming of the UvA
* Lesh, N., Mitzenmacher, M., & Whitesides, S. (2003). A complete and effective move set for simple protein folding. In Proceedings of the 7th annual international conference on research in computational molecular biology (RECOMB) (pp. 188–195). New York: ACM Press.
