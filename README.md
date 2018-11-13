# GrueLattice

## A Fully Homomophic Encryption Librairy

GrueLattice is an open-source software. The GrueLattice librairy is based on a Fully Homomorphic Encryption 
scheme and makes use of the [FFTW](http://www.fftw.org) librairy (the "Fastest Fourier Transform in the West").
GrueLattice name is based on the crane bird (grue in french), and on the term lattice. The librairy provides 
a symetric encryption scheme to encrypt and decrypt single bit messages.

The final version of GrueLattice is not planed to be, as for now, published.

** What is in the upper version **
- HMAC sinature for user authentication
- communication with the server
- lattice based blockchain to store and secured the data

### Requirements
GrueLattice requires the FFTW 3 librairy available at <http://www/fftw.org/download.html>, and a c++ compiler.
Also make sure that you have the valgrind library available at <https://github.com/smparkes/valgrind>.
The librairy is written primarily in C but a C++ compiler is needed to support a few synthical extensions (like
namespace) use to improve the readability of the code. 

### Check

Just run the command ``` make check```.

### Installation
To build the librairy, just run ```make```.

You can test the librairy by running the test program ```grueLattice```.

### Tests

In case you want to modify parameters and check correctness of the scheme, you can run the different unit tests in the `./tests` folder.

```
    cd tests/
    make
```

### Stats (please wait for this --> rewritting procedure in action)

In the `stats` folder, you have one program for each elementary operations and one global program that covers the whole scheme. These programs execute the operation(s) and measure the output error.
The input error variances for each operation is set, by default, to the expected output value of the operation that comes just before during a complete gate execution, for the given parameter sets.

```
    cd stats/
    make
    ./global
```

### Documentation

All the code source is commented with Doxygen.
In order to generate and access the documentation, just run ```sh doc.sh```. 
This will build the doxygen documentation and launch it in the firefox browser.
