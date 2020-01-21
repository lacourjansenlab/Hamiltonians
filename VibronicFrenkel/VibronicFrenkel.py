#
#   VibronicFrenkel.py
#
#       Description:
#           Generates a Frenkel exciton Hamiltonian with a Holstein like
#           vibrational mode coupled to it.
#
#       Usage:
#           python3 VibronicFrenkel.py [options]
#
#           for help type
#           python3 VibronicFrenkel.py --help
#
#       Dependencies:
#           Should work with any Python3 installation (for example Anaconda)
#               - uses: numpy, os, argparse, json, jsonschema, and unittest
#
#       Written by Marick Manrho
#       email: m.manrho@rug.nl
#
#       Last update: 15-jan-2020
#

#------------------------------------------------------------------------------#
# Hamiltonian generator
def gen_hamiltonian(p,ops,fc_table):
    import numpy as np

    # Initiate and generate
    H = np.zeros((p.size,p.size),dtype=np.double)
    H = gen_h_ops(p,fc_table,ops,H)

    if (p.v):
        print('Hamiltonian:\n', H)

    return H

#------------------------------------------------------------------------------#
# Index one-particle states
def index_ops(p):
    import numpy as np
    ops = np.zeros((p.N,p.MaxVib+1),dtype=np.int)-1
    count = 0
    for n in range(p.N):
        for v in range(p.MaxVib+1):
            ops[n,v] = count
            count = count + 1

    return ops

#------------------------------------------------------------------------------#
# Generate transition dipole moments
def gen_dipoles(p,ops,fc_table):
    import numpy as np
    mu = np.zeros((p.size,3))

    for n in range(p.N):
        for v in range(p.MaxVib+1):
            loca = ops[n,v]
            if (loca<0): continue
            mu[loca,:] = np.multiply(p.oscillator_strength,fc_table[0,v])

    if (p.v):
        print('Dipoles:\n', mu)

    return mu

#------------------------------------------------------------------------------#
# Generate gen_positions
def gen_positions(p,ops):
    import numpy as np
    pos = np.zeros((p.size,3))

    for n in range(p.N):
        for v in range(p.MaxVib+1):
            loca = ops[n,v]
            if (loca<0): continue
            pos[loca,0] = n

    if (p.v):
        print('Position vectors:\n', pos)

    return pos

#------------------------------------------------------------------------------#
# Generate one-particle states Hamiltonian
def gen_h_ops(p,fc_table,ops,H):
    import numpy as np

    # Diagonal elements
    for n in range(p.N):
        for v in range(p.MaxVib+1):
            loca = ops[n,v]
            if (loca < 0): continue
            randn = np.random.randn()
            H[loca,loca] = v*p.wvib + p.E + p.sigma*randn

    # Off diagonal elements
    for n in range(p.N):
        for v in range(p.MaxVib+1):
            loca = ops[n,v]
            if (loca < 0): continue
            for m in range(p.N):
                if (n == m): continue
                if (np.abs(n-m)>1): continue
                for q in range(p.MaxVib+1):
                    locb = ops[m,q]
                    if (locb < 0 or loca == locb): continue
                    H[loca,locb] = p.J*fc_table[v,0]*fc_table[0,q]
    return H

#------------------------------------------------------------------------------#
# Franck-Condon factor generator
def gen_fc_table(S,MaxVib):
    import numpy as np
    fc_table = np.zeros((MaxVib+1,MaxVib+1))
    for n in range(MaxVib+1):
        for m in range(MaxVib+1):
            fc_table[n,m] = gen_fc(n,m,S)

    return fc_table

def gen_fc(n,m,S):
    import numpy as np

    facn = factorial(n)
    facm = factorial(m)
    fc = 0
    for t in range(m+1):
        if (n-m+t < 0):
            pass
        else:
            facmt = factorial(m-t)
            facnmt = factorial(n-m+t)
            fact = factorial(t)
            facin = 1/(fact*facmt*facnmt)
            fc = fc + facin*S**(t*0.5)*S**(0.5*(n-m+t))*(-1)**(n-m+t)

    fc = fc*np.sqrt(facm*facn)*np.exp(-S/2.0)
    return fc

def factorial(n):
    f = 1
    for m in range(1,n+1):
        f = f*m
    return f

#------------------------------------------------------------------------------#
default_parameters = {
    "N": 10,
    "E": 0,
    "J": 100,
    "MaxVib": 5,
    "wvib": 1400,
    "S": 0.6,
    "sigma": 100,
    "oscillator_strength": [1,0,0],
    "units": "wavenumbers"
}

#------------------------------------------------------------------------------#
# Parameter parser
class Autovar(object):
    def __init__(self,data):
        self.__dict__ = data

#------------------------------------------------------------------------------#
# Engine which reads arguments, reads parameters, and run Hamiltonian generator
def main():
    # Look at arguments provided by user
    import argparse
    parser = argparse.ArgumentParser(prog='VibronicFrenkel', usage='python3 %(prog)s.py [options]')
    parser.add_argument('-I', default='default', help='Path to input file (json).')
    parser.add_argument('-v', action='store_true', help='Run verbose mode.')
    parser.add_argument('--test', action='store_true', help='Run test suite and exit.')
    parser.add_argument('--absorption', action='store_true', help='Plot absorption spectrum.')
    parser.add_argument('--structure', action='store_true', help='Visualize structure.')
    args = parser.parse_args()

    # Run tests
    if args.test:
        import os
        os.system('python3 Test_VibronicFrenkel.py')
        exit()

    # Load defaults or import json file
    if (args.I == 'default'):
        p = Autovar(default_parameters)
    else:
        import json
        import jsonschema

        input_file = open(args.I).read()
        import parameters_schema as schema

        # Validate input file against schema
        try:
            jsonschema.validate(json.loads(input_file), schema.schema)
        except jsonschema.ValidationError as e:
            print(e.message)
        except jsonschema.SchemaError as e:
            print(e)

        p = Autovar(json.loads(input_file))

    # Store verbose
    p.v = args.v

    # Determine number of states, for initiating arrays
    p.size = p.N*(p.MaxVib+1)

    # Generate Franck-Condon factors between ground state and excited state
    fc_table = gen_fc_table(p.S,p.MaxVib)

    # Index one-particle states
    ops = index_ops(p)

    # Generate dipole-moments
    mu = gen_dipoles(p,ops,fc_table)

    # Generate positions
    pos = gen_positions(p,ops)

    # Generate Hamiltonian
    H = gen_hamiltonian(p,ops,fc_table)

    # Plot absorption spectrum
    if (args.absorption):
        from Spectra.Spectra import absorption
        absorption(H,mu,p.size,p.sigma)

    # Plot structure
    if (args.structure):
        from Structure.Structure import visual
        visual(pos,mu,p.N,0.2)
        import matplotlib.pyplot as plt

#------------------------------------------------------------------------------#
# Run code
if __name__ == "__main__" and __package__ is None:
    # Ugly hack to allow absolute import from the root folder
    # whatever its name is. Please forgive the heresy.
    from sys import path
    from os.path import dirname as dir

    path.append(dir(path[0]))
    __package__ = "VibronicFrenkel"

    #path.append(dir(path[0]))
    #__package__ = "Structure"

if __name__ == '__main__':
    main()
