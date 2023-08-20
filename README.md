# OL-Tracker Command-Line Tool

To perform the analysis, clone this repo and run the following command inside the newly created folder in the terminal:

```
Rscript --vanilla origin-tracker.R [options]
```

Within the folder, there are two subdirectories namely ```/examples``` and ```/resources```.

The subdirectory ```/examples``` contains two files:

- ```arrival_dates.tsv``` this is TSV file where you put arrival time estimates, the required format should be self-explanatory
- ```query_countries.txt``` this is a simple TXT file where you can put countries that you are particularly interested in investigating as possible OLs

The subdirectory ```/resources``` contains two files:

- ```shortest_paths_Dmn_matrix_2020.renamed.csv``` this is CSV file containing the effective distances (extracted from shortest-path trees) with each country (with available air traffic data from 2020) as the presumed OL
- ```shortest_paths_Dmn_matrix_2021.renamed.csv``` equivalent to above but generated using air traffic data from 2021

## How-to

### Create ```arrival_dates.tsv```
It should be self-explainatory from the example how the estimated arrival times should be specified. To ensure consistency between country names, there is an option ```-s CHARACTER, --search=CHARACTER``` that allows you to check if the country you are adding has available air traffic data. Note that any specified countries without available air traffic data are ignored in the analysis. Make sure that you specify the path to file containing the estimated dates using the option ```-d CHARACTER, --dates=CHARACTER```.

Alternatively, you can also get a list of all countries with available data using the flag ```-n, --names```.

### Hypothesis testing - country X
Say you are interested in whether a specific country X is the OL. You are make it easier to assess this hypothesis by adding country X in the ```query_countries.txt``` file in the ```/examples``` folder. You can enter more than one countries here. Two additional files will be generated for these query countries to easier comparison and assessment. Make sure that you specify the path to file containing the list of query countries using the option ```-q CHARACTER, --query_countries=CHARACTER```.

### Full run

An example including options described above looks like

```Rscript --vanilla origin-tracker.R -d examples/arrival_dates.tsv -q examples/query_countries.txt```

#### Notes
- Run ```Rscript --vanilla origin-tracker.R -h``` to more details of available options