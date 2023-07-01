# PPSH-41
C++ for Computational Commutative Algebra(Like F2-Polynomial) and Finite Fields;

![logo](logo.gif)

How to use:

1. Compiling the source code;
```bash
$ g++ ppsh_functions.cpp -lginac -lcln -fPIC -shared -o libppsh.so
```

2. Compiling the your own code (here's the testing code):
```bash
$ g++ test_ppsh_functions.cpp -lginac -lcln -L. -lppsh -o exe
$ ./exe
g^1 = 1+a+2*a^2
g^2 = 2*a+a^2
g^3 = 1+2*a+a^2
g^4 = 2
g^5 = 2+2*a+4*a^2
g^6 = a+2*a^2
g^7 = 2+a+2*a^2
g^8 = 1
g^9 = 1+a+2*a^2
```
