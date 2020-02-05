# beerprot

Part one:
  
  1.1 BLAST <br>
        1.1.a  From starting sequence we used BLAST for retrieving similar sequences <br>
        1.2.b  Modified the BLAST results taking out outliers <br>
  1.2. Multiple Sequence Alligment (MSA) <br>
        1.2.a  Generated a MSA from the BLAST results using CLUSTELW(see 1.1.a) <br>
        1.2.b  Modified MSA with Jelview <br>
  1.3. Hidden Markov Model (HMM) <br>
        1.3.a  Build an HMM model using the MSA of point 1.2.b <br>
        1.3.b  Retireve a set of human proteins with the same PFAM of the initial sequence (positive examples) <br>
        1.3.c  Retrieve a set of human proteins with a different PFAM of the initial sequence (negative examples) <br>
        1.3.d  Evaluated the HMM model using the positive and the negative examples <br>
