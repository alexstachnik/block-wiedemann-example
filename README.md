block-wiedemann-example
=======================

A Python implementation of the Block-Wiedemann (Coppersmith) algorithm.

Warning: This code is buggy and painfully slow. It was never meant for
serious use, just for reference.

Look at main.py for a usage example. It generates a non-singular
matrix and computes the minmial polynomial by Wiedemann's algorithm
(Solving sparse linear equations over finite fields, D. H. Wiedemann,
IEEE Trans. Inform. Theory (1986)) and by Coppersmith's
Block-Wiedemann algorithm (Solving Homogeneous Linear Equations Over
GF(2) Via Block Wiedemann Algorithm, Don Coppersmith, Mathematics of
Computations (1994)).

