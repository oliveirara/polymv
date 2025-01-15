# polymv_C.pyx
from libc.stdlib cimport malloc, free
cimport numpy as np
import numpy as np

cdef extern from "mvs.h":
    void multipol_vec(double* input1, double* input2, double* output1, double* output2, int LMAX)

cdef extern from "fvs.h":
    void frechet_vec(double* input1, double* input2, double* output1, double* output2, int LMAX)

def mvs(np.ndarray[double, ndim=1] input1, np.ndarray[double, ndim=1] input2, int LMAX):
    
    # Garante que os arrays são contíguos
    input1 = np.ascontiguousarray(input1, dtype=np.float64)
    input2 = np.ascontiguousarray(input2, dtype=np.float64)

    cdef int tamanho_output = ((LMAX * (LMAX + 1)) - 2)
    
    cdef np.ndarray[double, ndim=1] output1 = np.zeros(tamanho_output, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] output2 = np.zeros(tamanho_output, dtype=np.float64)


    multipol_vec(&input1[0], &input2[0], &output1[0], &output2[0], LMAX)

    return output1, output2


def fvs(np.ndarray[double, ndim=1] input1, np.ndarray[double, ndim=1] input2, int LMAX):
    
    # Garante que os arrays são contíguos
    input1 = np.ascontiguousarray(input1, dtype=np.float64)
    input2 = np.ascontiguousarray(input2, dtype=np.float64)

    cdef int tamanho_output = LMAX - 2 
    
    cdef np.ndarray[double, ndim=1] output1 = np.zeros(tamanho_output, dtype=np.float64)
    cdef np.ndarray[double, ndim=1] output2 = np.zeros(tamanho_output, dtype=np.float64)

  
    frechet_vec(&input1[0], &input2[0], &output1[0], &output2[0], LMAX)

    return output1, output2



