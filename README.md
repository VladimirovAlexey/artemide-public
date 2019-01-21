# artemide-public
The public repository of artemide package for TMD-physics (transverse momentum dependent)

------------------------------------------------------------------------------------------------------
    CHECK:
    In makefile set (in the begining of file)
    FCompiler        <= your prefered fortran compiler (f95 at least, gfortran also works)
    Fflags        <= if you use openmp
    FOPT        <= For extra options, links,etc. see LHAPDF
    
    the harpy compiles with the help of f2py package from numpy (python2)
    
    The file "constants" must be in the same location as your program. Check it, it accumulates all options.


------------------------------------------------------------------------------------------------------
    LHAPDF:
    By default, artemida uses LHAPDF (although you can put your own PDF's, see manual).
    So, make sure that LHAPDF properly installed and check the link to it.
    For artemide)
        in FOPT of makefile link to LHAPDF: it should look like
        FOPT = -L/path/to/LHAPDF/Installation/lib -lLHAPDF -lstdc++
    
    For harpy)
        make sure that PYTHONPATH contains path for LHAPDF it typically looks like
        /path/to/LHAPDF/Installation/lib/python2.7/site-packages

------------------------------------------------------------------------------------------------------
Commands in make

make
=> Compiles the artemide package (necessary for all others)

make harpy
=> Compiles the harpy from artemide

make test
=> Compiles and execute test.f90. It initializes artemide and calculates some cross-section. Also tests parallel computation (if it is)

make program TARGET=abc.f90
=> Compiles a program abc.f90 in the /Prog with artemide


-------------------------------------------------------------------------------------------------------
See manual for details on artemide  in /doc
If you have questions, suggestions => E-mail: vladimirov.aleksey@gmail.com


