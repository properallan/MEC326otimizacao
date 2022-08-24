# Intallation

Clone this repository and goes to the it's root path. To compile the eulerQ1D solver using gcc, just type on a linux terminal 

```./compile.sh eulerQ1D.cpp```

** Note: `compile.sh` must have execution permision, needed only at first compile (`chmod +x compile.sh`)

## Requirements

All dependencies are installable by

```pip install -r requirements.txt```

# Usage

Sample line of code for an iverse design using the polynomial parametrization (runs a lot faster than with splines)

```python optXti.py -p 1 -m trust-constr -wd ./resultsTrustBounds/ -r ./resTrustBounds.pickle -tol 1e-4```
