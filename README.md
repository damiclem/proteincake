# beerprot

## Part one: Models
This part explains every model developed in this project and how to use them. Some models, such as PSSM or HMMER require an initial MSA, which can be retrieved by following the instructions above, while others, such as JACKHMMER, don't.

*Note* that each model described here returns a table with two columns:
1. Protein accession number;
2. Set of amino-acid positions in which the domain has been found.

### 1) Position Specific Scoring Matrix (PSSM)
This model is actually a PSI-BLAST algorithm which takes as input a specific PSSM. Parameters are:
- *fit*: defines wether to fit or not a new model. Boolean, default *True*;
- *blast_path*: path to blast file, if fit is True. String, default *data/blast.fasta*;
- *msa_path*: path to multiple sequence alignment result, FASTA formatted. String, default *msa.edited.fasta*;
- *model_path*: path where to store the new model, if fit is True. String, default *models/model.pssm*
- *test_path*: path from which test dataset is loaded. String, default *data/human.fasta*;
- *out_path*: path where to store model results. String;
- *num_iterations*: number of PSI-BLAST iterations. Int, default *3*;
- *e_vale*: e-value threshold on results. Float, default *0.001*.

```shell
python modules/pssm.py --fit True --blast_path path/to/blast.fasta --msa_path path/to/msa.fasta --model_path path/to/model.pssm --test_path path/to/test.fasta --out_path path/to/out.tsv --num_iterations 3 --e_value 0.001
```

### 2) Hidden Markov Model
This model allows to use either HMMER or JACKHMMER in order to obtain domain classification and domain's positions matching. Parameters are:
- *algorithm*: whether to use HMMER or JACKHMMER. String, default *hmmer*;
- *fit*: wheter to fit a new model or not. Bool, default *True*;
- *seq_path*: path to query sequence, in case selected algorithm is JACKHMMER. String, default *data/domain.fasta*;
- *msa_path*: path to multiple sequence alignment, in case selected algorithm is HMMER. String, default *data/msa.edited.fasta*;
- *test_path*: path to FASTA formatted test set. String, default *data/human.fasta*;
- *model_path*: path to model to load/create, in case selected algorithm is HMMER. String, default *models/model.hmm*;
- *out_path*: path where to store model output. String;
- *e_value*: e-value threshold. Float, default *0.001*;

```shell
python modules/hmm.py --algorithm hmmer --fit True --seq_path path/to/query/sequence.fasta --msa_path path/to/msa.fasta --test_path path/to/test.fasta --model_path path/to/model --out_path path/to/out.tsv --e_value
```

### 3) Ensemble model
Differently from PSSM and HMM models, ensemble model takes as input the output of other models and runs majority voting after clusterizing predicted domains. Parameters are:
- *models_out*: list of models outputs. List;
- *out_path*: path where to store the model result. String.


```shell
python modules/ensemble.py --models_out path/to/model1.tsv path/to/model2.tsv --out_path path/to/out.tsv
```

## Part Two: Annotation Enrichment
In this section we have performed Enrichment, using first the **GO** and secondly the **DO** terms, on all the four kind of datasets: *Original*, *PDB Network*, *Architectures*, *STRING Network*. All the results are inside the paths *results/go_enrichment* and *results/do_enrichment*.

### Method:
The enrichments, measured by using the *Fisher Test*, is done by measuring the **Odd Ratio** of the presence of a *GO/DO* w.r.t two datasets, in other words by comparing how much a *GO/DO* term is present in a dataset, the *target dataset*, w.r.t. another dataset, the *background dataset*.

An Odd Ratio close to 0 indicates that a *GO/DO* term is more present in the *background dataset*, a value close to 1 indicates that the *GO/DO* term is present in equal measure in between the two datasets and a value significantly higer than 1 (potentially *plus infinity*) indicates that the term is present more in the *target dataset* (please note that *a term is more present* in a dataset w.r.t. another dataset means that the ratio of that term in a dataset is higher than the ratio of the other dataset, not that a term has been found more times in a dataset w.r.t. another one).

After computing the *OddRatio* and the *p-value* of every *GO/DO* term, which was present in both the *target* and *background* dataset, we computed the depth of this terms using the *Gene/Diesease Ontology Graph*, we then take out all the terms with a too high *p-value* or with a too high depth (these thresholds are parameters that can differ from dataset to dataset, deafult parameters are: <br>
- *maximum p-value*: 0.05 <br>
- *maximum depth*: 5

Moreover, we decided to transmit the p-value of a children to its parents, in order, to easier the presence of *high level terms* in the final results.

For each *GO/DO* term that passed the filter, we annoted its *GO/DO id*, *p-value*, *depth*, *description* and *Score*, with Score being the natural logarithm of 1 over p-value, *ln(1/p-value)*. The score is used to generate the *WordCloud*, such that the description of the GO/DO terms with lower *p-value* would appear with a bigger font.  

### 1) Original
The enrichment has been done by using the following datasets: <br>
- *Target*: Original <br>
- *Background*: All human proteome on SwissProt
 
 The code used for displaying the GO result is:
 
    - python modules/go_modules/enrichment_go.py
    
### 2) Architectures
For the enrichment this time we have a set of *target datasets*, one for every architecture, where an architecture is characterized by the unordered set of *PFAM domains*. In summary: <br>
- *Target*: Set of original proteins that share the same PFAM domains (e.g. PF00169, PF00397) <br>
- *Background*: Original dataset.

The code used for creating all the architecture datasets and performing the enrichment on one single target dataset (in this case the dataset characterized by all the proteins belonging only to *PF00397*) is:

    - python modules/go_modules/architecture.py
    - python modules/go_modules/enrichment_go.py  --target_path "data/architecture/go_architectures/PF00397_arch.csv" --background_path "data/architecture/go_architectures/architecture_background.csv"
    
If you want to run the enrichment tests on all the *target datasets* you can run:

    - python modules/go_modules/architecture.py
    - python modules/go_modules/test_all_architectures.py
    
Be aware that this script will not display the results but it will save them in a folder (*default: results/go_enrichment/architectures*). Please see the scipt to see all the arguments.
    

### 3) PDB Network
The enrichment has been done by using the following datasets: <br>
- *Target*: Original (with a PDB) plus other human proteins which are found as chain in the same PDB <br>
- *Background*: All human proteome on SwissProt with a PDB
 
 The code used for creating the PDB datasets and displaying the enrichment results is:
 
    - python modules/go_modules/pdb_network.py
    - python modules/go_modules/enrichment_go.py  --target_path "data/pdb/pdb_target_go.csv" --background_path "data/pdb/pdb_background_go.csv"

### 4) STRING Network
For the STRING part, the enrichment has been done by using the following datasets: <br>
- *Target*: Original plus the direct interactors found in the STRING database <br>
- *Background*: All human proteome on SwissProt intersected with the STRING database.

Note that we considered a protein as a direct interactor with one in the original dataset if the **Combined Score** is higher than 700, in other words, only if the probability of observing an interaction between these two proteins, taking into consideration the *prior probability* of observing a random interaction, is above 0.70.    
 
 The code used for creating the STRING datasets and displaying the enrichment results is:
 
    - python modules/go_modules/string_network.py
    - python modules/go_modules/enrichment_go.py  --target_path "data/string/string_target_go.csv" --background_path "data/string/string_background_go.csv" 
    
### Note on the Enrichment part
Note 1: <br>
The code showed above produced is an example of the *GO* part, to perform the *DO* part just change the *go* to *do*. E.g:

    - python modules/do_modules/string_network.py
    - python modules/do_modules/enrichment_do.py  --target_path "data/string/string_target_do.csv" --background_path "data/string/string_background_do.csv"

Note 2: <br>
There are two *jupyter notebook*, *GO_demo.ipynb* and *DO_demo.ipynb* that show how all the dataset are generated and how the results are obtained.

Note 3: <br>
Additional arguments, that are not shown in the code, are present in the modules, please see the code to see all the arguments. <br>
Example of arguments: <br>
- *ontology_path*: the path of the *GO/DO* ontologies graph files. Default for *GO* is: data/go/go.json.gz; <br>
- *out_path*: the path where to save the results (*pandas.dataframe*). Deafualt is None (print the output, do not save it); <br>  
- *out_wordcloud*: the path where to save the *WorldCloud* (*matplotlib.pyplot.figure*). Deafualt is None (show the output, do not save it); <br>
- *p_value*: maximum *p-value* that a *GO/DO* term must have in order to be considered significant. Dedault is 0.05; <br>
- *depth*: maximum *depth* on the ontology graph that a *GO/DO* term must have in order to be considered significant. Dedault is 4; <br>
- *bonferroni*: Perform *Bonferroni correction* (1) or dont't (0). Default is 1. 
