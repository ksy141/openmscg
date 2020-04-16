:orphan:

CGTOP File Format
=================

Provide a Gromacs-style topology file for the CG system with filename 'top.in'. The
bonding described in this file is the bonding that will be used to set up all bonded
interaction potential basis sets. Bond types are automatically inferred from site types
and cannot be specified any other way.

Example inputs can be found in the subdirectory of the example directory. The general
format of the topology is outlined below::

    cgsites Ncg (* Ncg =Total number of coarse-grained sites)
    cgtypes Nt  (* Nt = Total number of coarse-grained site types)
    ( Site type names, row number corresponding to type number)
    ...
    ...
    (Additional add-ons go here for density interactions)
    moltypes Nm (* Nm = Total number of molecule types)
    mol Ns b    (Header for the block defining the first molecule type)
                (* Ns = number of coarse-grained sites in this molecule)
                (* b = bonded interaction specification style)
                (** b determines if angles and dihedrals are to be input manually or inferred
                    from the pair bond topology)
                (** b = 3 means bonds/angles/dihedrals are inferred from the pair bond list
                    i.e., 1-2,1-3,1-4 interactions are not included in the force matching of
                    nonbonded interactions)
                (** b = 2 means only angles are inferred)
                (** b = 1 means no interactions beyond pair bonds are present in the molecule)
                (** b = -1 means angles and dihedrals are manually specified in this file:
                    the sections "angles" and "dihedrals" need to be provided)
    sitetypes   (Define a site type for each site in the molecule)
    t           (* t = type of this site)
                (Note: only 1 type per line; type numbers are 1-based)
    ...         (Should have a type for all Ns sites in this molecule)
    bonds Nb    (Header for the block defining all pair bonds between sites in the molecule)
                (* Nb  = number of pair bonds)
    i j         (For each bond i-j, site indices i and j use the ordering in sitetypes)
    . .         (Should have entries for all Nb bonds in this molecule)
    angles  Na Fa (Header for the block defining manually-specified angles between sites in
                 the molecule)
                (* Na = number of angles)
                (* Fa = angle format)
    i j k		( ** Fa = 1 : For each angle i-j-k, site indices i, j, and k use the ordering in sitetypes)
    j i k       ( ** Fa = 0 : For each angle i-j-k, site indices i, j, and k use the ordering in sitetypes)
                ( Note: the order of sites with Fa = 0 specifies the central site j first)
    . . .       (Should have entries for all Na angles in this molecule)
    dihedrals Nd Fd (Header for the block defining the manually-specified dihedrals between sites
                  in the molecule)
                ( * Nd = number of dihedrals
                ( * Fd = dihedral format)
    i j k l		( ** Fd = 1 : For each dihedral i-j-k-l, site indices i, j, k, and l use the ordering in sitetypes)
    j k i l     ( ** Fd = 0 : For each dihedral i-j-k-l, site indices i, j, k, and l use the ordering in sitetypes)
                ( Note: the order of sites with Fd = 0 specifies the central sites j and k first)
    . . . .     (Should have entries for all Nd dihedrals in this molecule)
    . . . . . . (Define additional molecule types until Nm molecule types have been defined)
                (Note: this goes from the mol tag through the bonds, angles, or dihedral tags)
    system s    (Header for the system definition block)
                (* s = states number of groups of molecules comprising the system)
    Mt Mn       (* Mt = The molecule type using the ordering of the mol labels above)
                (* Mn = The number molecules of this type that continuous in frame ordering)
    ... ...     (Should have s rows of data)