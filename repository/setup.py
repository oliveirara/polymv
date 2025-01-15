from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

exts = [
    Extension(
        "polymv_C",
        sources=["polymv_C.pyx","mvs.c", "fvs.c"],
        libraries=["m", "gmp", "mps",  "chealpix", "cfitsio", "nlopt"],             # Adiciona a biblioteca GMP
#        library_dirs=["/usr/lib", "/usr/local/src/cfitsio-4.4.0", "/usr/local/src/Healpix_3.82/lib", ],      # Diretório onde a GMP está instalada
        library_dirs=["/usr/lib", "/usr/local/lib"],      # Diretório onde a GMP está instalada
        include_dirs=[np.get_include()],
#        include_dirs=[np.get_include(), "/usr/local/src/cfitsio-4.4.0", "/usr/local/src/Healpix_3.82/include", "/usr/local/include/nlopt"], # Adiciona os diretórios de inclusão
        extra_compile_args=["-std=c99", "-fopenmp", "-fPIC", "-march=native", "-ffast-math", "-mfpmath=sse", "-O2"], # Compatibilidade com C99 e OpenMP
        extra_link_args=['-fopenmp', '-lgomp'],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )
]

setup(
    ext_modules=cythonize(exts),
)