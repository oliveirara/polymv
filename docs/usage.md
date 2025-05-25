# polymv Usage Guide

This guide demonstrates how to use the main features of the `polymv` package for working with multipole vectors.

---

## Setup

Import the necessary libraries:

```python
import numpy as np
import healpy as hp

from polymv import polymv
from polymv import utils
```

---

## Generating Test Data

Let's create a test set of $a_{\ell m}$s:

```python
lmax = 1500
np.random.seed(1)
alm = hp.synalm(np.ones(2000), lmax=lmax)
alm = np.column_stack([alm.real, alm.imag])
alm.shape
```

**Output:**
```
(1127251, 2)
```

---

## Multipole Vectors (`mvs`)

To obtain multipole vectors from $\ell=2$ to $\ell=1500$:

```python
mvs_1500 = polymv.mvs(alm[:, 0], alm[:, 1], lmax)
np.vstack(mvs_1500).shape
```

**Output:**
```
(2, 2251498)
```

You can view a sample of the multipole vectors:

```python
np.vstack(mvs_1500)[:, :5]
```

**Output:**
```
array([[ 1.83441666,  1.37302603,  1.30717599,  1.52087956,  1.62071309],
       [ 0.12767786,  1.33828394, -3.01391479,  0.19300013, -2.94859252]])
```

---

## Utility Functions

### Extracting Specific Multipole Vectors

Get multipole vectors for $\ell=3$:

```python
utils.get_mvs_from_many(np.vstack(mvs_1500).T, 2, lmax, 3)
```

**Output:**
```
array([[ 2.73836992,  2.08191869],
       [ 1.52087956,  0.19300013],
       [ 1.62071309, -2.94859252],
       [ 2.34022757,  0.13144151],
       [ 0.40322273, -1.05967396],
       [ 0.80136509, -3.01015115]])
```

---

### Convert to Cartesian Coordinates

Convert the multipole vectors to Cartesian coordinates:

```python
cart_coords = utils.to_cart(np.vstack(mvs_1500).T)
cart_coords[:3]
```

**Output:**
```
array([[ 0.95759438,  0.12293233, -0.26057752],
       [ 0.22593142,  0.95412221,  0.19648358],
       [-0.95759438, -0.12293233,  0.26057752]])
```

---

### Filter for North Hemisphere

Get only vectors in the north hemisphere:

```python
north_mvs = utils.mvs_north(np.vstack(mvs_1500).T)
north_mvs.shape
```

**Output:**
```
(1125749, 2)
```

---

## Fréchet Mean Vectors (`fvs`)

Compute Fréchet mean vectors for all multipoles:

```python
fvs = polymv.fvs(mvs_1500[0], mvs_1500[1], lmax)
np.vstack(fvs).shape
```

**Output:**
```
(2, 1498)
```

Show a sample:

```python
np.vstack(fvs)[:, :5]
```

**Output:**
```
array([[0.39372199, 1.30678037, 1.41442982, 1.18631234, 1.2620573 ],
       [5.54727987, 1.64990686, 0.53026658, 2.88028155, 0.76763323]])
```
