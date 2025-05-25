# polymv Test Example

This example demonstrates the usage and performance of the `polymv` package on a DELL XPS 8900 with a 6th Generation Intel® Core™ i5-6400 Processor (2.7 GHz up to 3.3 GHz, 6MB Cache).

---

## Setup

Import the necessary libraries:

```python
import numpy as np
import healpy as hp

from polymv import polymv
```

---

## Explore the `polymv` Module

Check the available functions and their docstrings:

```python
help(polymv)
```

**Output:**
```
Help on module polymv.polymv in polymv:

NAME
    polymv.polymv

FUNCTIONS
    fvs(input1, input2, LMAX)
        Get Fréchet mean vectors in spherical coordinates for a multipole given a map based on lmax.

        Parameters:
            input1 (float): Multipole vector theta coordinate.
            input2 (float): Multipole vector phi coordinate.
            LMAX (int): Maximum l value to compute the Fréchet mean.

        Returns:
            array: An array for all Fréchet mean based on all multipole vectors in spherical coordinates (theta, phi) from l=2 to lmax.

    mvs(input1, input2, LMAX)
        Get multipole vectors in spherical coordinates for a multipole given a map based on lmax.

        Parameters:
            input1 (float): Real part of multipole moments from a map (Healpy indexing).
            input2 (float): Imaginary part of multipole moments from a map (Healpy indexing).
            LMAX (int): Maximum l value to compute the multipole vectors.

        Returns:
            array: An array for all multipole vectors in spherical coordinates (theta, phi in radians) from l=2 to lmax.

DATA
    __test__ = {}

FILE
    /home/renan/polymv-project/polymv/polymv/polymv.cpython-312-x86_64-linux-gnu.so
```

---

## Set Maximum Multipole

```python
lmax = 1000
```

---

## Generate Random Multipole Data

```python
%%time
np.random.seed(1)
alm = hp.synalm(np.ones(2000), lmax=lmax)
alm = np.column_stack([alm.real, alm.imag])
alm.shape
```

**Output:**
```
CPU times: user 43 ms, sys: 11.8 ms, total: 54.8 ms
Wall time: 54.3 ms

(501501, 2)
```

---

## Compute Multipole Vectors

```python
%%time
mvs = polymv.mvs(alm[:, 0], alm[:, 1], lmax)
print(type(mvs), [arr.shape for arr in mvs])
```

**Output:**
```
CPU times: user 5min 50s, sys: 3.91 s, total: 5min 53s
Wall time: 1min 37s

<class 'tuple'> [(1000998,), (1000998,)]
```

Show a sample of the multipole vectors:

```python
print("theta:", mvs[0][:5])
print("phi:", mvs[1][:5])
```

**Output:**
```
theta: [1.44757865 1.13235582 1.694014   1.07304038 0.03596525]
phi:   [-0.23422314  1.99653803  2.90736952 -0.67228949  0.58064643]
```

---

## Compute Fréchet Mean Vectors

```python
%%time
fvs = polymv.fvs(mvs[0], mvs[1], lmax)
print(type(fvs), [arr.shape for arr in fvs])
```

**Output:**
```
CPU times: user 2min 46s, sys: 296 ms, total: 2min 47s
Wall time: 1min 10s

<class 'tuple'> [(999,), (999,)]
```

Show a sample of the Fréchet mean vectors:

```python
print("theta:", fvs[0][:5])
print("phi:", fvs[1][:5])
```

**Output:**
```
theta: [0.61105855 1.18990493 1.3201296  1.56128225 0.35433486]
phi:   [4.30053475 2.58923118 4.61933201 6.0707641  0.62762944]
```
