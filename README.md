# beerprot

## Part one:
  
  1.1 BLAST
  
        1.1.a  From starting sequence we used BLAST for retrieving similar sequences
        1.2.b  Modified the BLAST results taking out outliers
  
  1.2. Multiple Sequence Alligment (MSA)
  
        1.2.a  Generated a MSA from the BLAST results using CLUSTELW(see 1.1.a)
        1.2.b  Modified MSA with Jelview
  
  1.3. Hidden Markov Model (HMM)
        
        1.3.a  Build an HMM model using the MSA of point 1.2.b
        hmmbuild models/model.hmm data/msa.fasta
        1.3.b  Retrieve a set of human proteins with the same PFAM of the initial sequence (positive examples)
        database:(type:pfam pf00397) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"
        1.3.c  Retrieve a set of human proteins with a different PFAM of the initial sequence (negative examples)
        reviewed:yes AND organism:"Homo sapiens (Human) [9606]" AND NOT database:(type:pfam pf00397)
        1.3.d  Evaluated the HMM model using the positive and the negative examples
        hmmsearch --domtblout results/HMM_domtblout_pos -o results/HMM_results_pos models/model.hmm data/positive_pfam.fasta


## Part Two

### 2.1 Annotation Enrichment
In this section we have performed Enrichment, using first the **GO** and secondly the **DO** terms, on all the four kind of datasets: *Original*, *PDB Network*, *Architectures*, *STRING Network*.

**Method:** <br>
The enrichments, measured by using the *Fisher Test*, is done by measuring the **Odd Ratio** of the presence of a *GO/DO* w.r.t two datasets, in other words by comparing how much a *GO/DO* term is present in a dataset, the *target dataset*, w.r.t. another dataset, the *background dataset*.

An Odd Ratio close to 0 indicates that a *GO/DO* term is more present in the *background dataset*, a value close to 1 indicates that the *GO/DO* term is present in equal measure in between the two datasets and a value significantly higer than 1 (potentially *plus infinity*) indicates that the term is present more in the *target dataset* (please note that *a term is more present* in a dataset w.r.t. another dataset means that the ratio of that term in a dataset is higher than the ratio of the other dataset, not that a term has been found more times in a dataset w.r.t. another one).

After computing the *OddRatio* and the *p-value* of every *GO/DO* term, which was present in both the *target* and *background* dataset, we computed the depth of this terms using the *Gene/Diesease Ontology Graph*, we then take out all the terms with a too high *p-value* or with a too high depth (these thresholds are parameters that can differ from dataset to dataset, deafult parameters are: <br>
- *maximum p-value*: 0.05 <br>
- *maximum depth*: 5

Moreover, we decided to transmit the p-value of a children to its parents, in order, to easier the presence of *high level terms* in the final results.

For each *GO/DO* term that passed the filter, we annoted its *GO/DO id*, *p-value*, *depth*, *description* and *Score*, with Score being the natural logarithm of 1 over p-value, *ln(1/p-value)*. The score is used to generate the *WordCloud*, such that the description of the GO/DO terms with lower *p-value* would appear with a bigger font.  

#### 2.1.a Original
The enrichment has been done by using the following datasets: <br>
- *Target*: Original <br>
- *Background*: All human proteome on SwissProt
 
 The code used for displaying the GO result is:
 
    - python modules/go_modules/enrichment_go.py
    
#### 2.1.b Architectures
For the enrichment this time we have a set of *target dataset*, one for every architecture, where an architecture is characterized by the unordered set of *PFAM domains*. In summary: <br>
- *Target*: Set of original proteins that share the same PFAM domains (e.g. PF00169, PF00397) <br>
- *Background*: Original dataset.

The code used for creating all the architecture datasets and performing the enrichment on one single target datsset (in this case the dataset characterized by all the proteins belonging only to *PF00397*) is:

    - python modules/go_modules/architecture.py
    - python modules/go_modules/enrichment_go.py  --target_path "data/architecture/go_architectures/PF00397_arch.csv" --background_path "data/architecture/go_architectures/architecture_background.csv"    

#### 2.1.c PDB Network
The enrichment has been done by using the following datasets: <br>
- *Target*: Original (with a PDB) plus other human proteins which are found as chain in the same PDB <br>
- *Background*: All human proteome on SwissProt with a PDB
 
 The code used for creating the PDB datasets and displaying the enrichment results is:
 
    - python modules/go_modules/pdb_network.py
    - python modules/go_modules/enrichment_go.py  --target_path "data/pdb_data/pdb_target_go.csv" --background_path "data/pdb_data/pdb_background_go.csv"

#### 2.1.d STRING Network
For the STRING part, the enrichment has been done by using the following datasets: <br>
- *Target*: Original plus the direct interactors found in the STRING database <br>
- *Background*: All human proteome on SwissProt intersected with the STRING database.

Note that we considered a protein as a direct interactor with one in the original dataset if the **Combined Score** is higher than 700, in other words, only if the probability of observing an interaction between these two proteins, taking into consideration the *prior probability* of observing a random interaction, is above 0.70.    
 
 The code used for creating the STRING datasets and displaying the enrichment results is:
    - python modules/go_modules/string_network.py
    - python modules/go_modules/enrichment_go.py  --target_path "data/string/string_target_go.csv" --background_path "data/string/string_background_go.csv" 
    
#### 2.1.e: Note on the Enrichment part
Note 1: <br>
The code showed above produced is an example of the *GO* part, to perform the *DO* part just change the *go* to *do*. E.g:

    - python modules/do_modules/string_network.py
    - python modules/do_modules/enrichment_do.py  --target_path "data/string/string_target_do.csv" --background_path "data/string/string_background_do.csv"

Note 2: <br>
There are two *jupyter notebook*, *GO_demo.ipynb* and *DO_demo.ipynb* that show how all the dataset are generated and how the results are obtained.

Note 3: <br>
Additional arguments, that are not shown in the code, are present in the modules, please see the code to see all the arguments. <br>
Example of arguments: <br>
- ontology_path: the path of the *GO/DO* ontologies graph files. Default for *GO* is: data/go/go.json.gz; <br>
- out_path: the path where to save the results (*pandas.dataframe*). Deafualt is None (print the output, do not save it); <br>  
- out_wordcloud: the path where to save the *WorldCloud* (*matplotlib.pyplot.figure*). Deafualt is None (show the output, do not save it); <br>
- p_value: maximum *p-value* that a *GO/DO* term must have in order to be considered significant. Dedault is 0.05; <br>
- depth: maximum *depth* on the ontology graph that a *GO/DO* term must have in order to be considered significant. Dedault is 4; <br>

