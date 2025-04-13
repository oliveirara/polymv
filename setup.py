import os

import numpy as np
from setuptools import Extension, setup

home_dir = os.getenv("HOME", "")

setup(
    name="polymv",
    ext_modules=[
        Extension(
            "polymv.polymv",
            sources=[
                "src/polymv.pyx",
                "src/mvs/mvs.c",
                "src/fvs/fvs.c",
            ],
            libraries=[
                "gmp",
                "mps",
                "m",
                "cfitsio",
                "chealpix",
                "nlopt",
                "gomp",
            ],
            library_dirs=[
                os.path.join(home_dir, ".local", "lib"),
                os.path.join(home_dir, ".local", "lib64"),
            ],
            include_dirs=[
                np.get_include(),
                os.path.join(home_dir, ".local", "include"),
            ],
            extra_compile_args=[
                "-fPIC",
                "-Wall",
                "-march=native",
                "-O3",
                "-fopenmp",
                "-ffast-math",
                "-mfpmath=sse",
            ],
            extra_link_args=[
                "-lgmp",
                "-lgmpxx",
                "-lmps",
                "-lm",
                "-lcfitsio",
                "-lchealpix",
                "-lstdc++",
                "-lnlopt",
            ],
            define_macros=[
                ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ],
        ),
        Extension(
            "polymv.utils",
            sources=["src/utils.pyx"],
            include_dirs=[np.get_include()],
        ),
    ],
)
