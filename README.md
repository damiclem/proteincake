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

**1. Perform Enrichment on GO terms** <br>
Here is an example of the *modules* to launch to create all the dataset and to perform enrichment on them. <br>
To compare the 4 couple of datasets we first generate the data (in the first point the data is already generated and the path of it is set as deafult arguments), then we run the enrichment.

#### First point
>python modules/go_modules/enrichment_go.py

#### Second point
>python modules/go_modules/pdb_network.py <br>
>python modules/go_modules/enrichment_go.py  --target_path "data/pdb_data/pdb_target_go.csv" --background_path "data/pdb_data/pdb_background_go.csv"

#### Third point
>python modules/go_modules/architecture.py <br>
>python modules/go_modules/enrichment_go.py  --target_path "data/architecture/go_architectures/PF00397_arch.csv" --background_path "data/architecture/go_architectures/architecture_background.csv"

#### Fourth point
>python modules/go_modules/string_network.py <br>
>python modules/go_modules/enrichment_go.py  --target_path "data/string/string_target_go.csv" --background_path "data/string/string_background_go.csv"
