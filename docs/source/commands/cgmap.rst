cgmap
=====

.. include:: common.rst

.. automodule:: mscg.cli.cgmap

* |cli-opt-top-name|

* |cli-opt-traj|


CG Mapping File
---------------

The CG mapping file defines the rules of constructing CG coordinates from
all-atom (AA) trajectories. The file should be written in `YAML
<https://yaml.org/start.html>`__ format. An example of the CG mapping file is
below:

.. code-block::
   :caption: Table 1: CG mapping of the water box
   :emphasize-lines: 1,6
   
   site-types:
     WAT:
       index:    [   0,   1,   2]
       x-weight: [16.0, 1.0, 1.0]
       f-weight: [ 1.0, 1.0, 1.0]
   system:
     - anchor: 0
       repeat: 256
       offset: 3
       sites:
         - [WAT, 0]


This example demonstrates the CG mapping from a three-body water box to an
one-site topology. We assume that the AA configuration has 768 atoms as follows:

.. code-block::
   :caption: Table 2: AA configuration of the water box
   
       ID  Type  X   Y   Z
   1   OW    ... ... ...
   2   HW    ... ... ...
   3   HW    ... ... ...
   4   OW    ... ... ...
   5   HW    ... ... ...
   6   HW    ... ... ...
   ... ..    ... ... ...

    
The content in the mapping YAML file is parsed by Python as key-value pairs
(like a python dictionary). The CG mapping file should have two top-level
sections: `site-types` and `system`, which are described below.

Site Types
^^^^^^^^^^

Each entry in the `site-types` section defines mapping from multiple atoms to a single CG site labeled by the name of the entry. In the example above, there's only one CG type (`WAT`), each site of which will be mapped from three atoms.

Each CG site entry has three sub-entries, `index`, `x-weight`, and `f-weight`, which must be in the same size as the number of atoms mapping to the CG site.

* The `index` list defines the index atoms to be used. The index are the `offsets` to an **anchor** atom, which will be described later. In this example, three atoms in a sequence (0, 1 and 2) relative to the anchor atom will be used for mapping.

* The `x-weight` list defines the weighting factors to map the CG **coordinates** by the equation below, in which the subscript *i* goes over all involved atoms for this CG site. Please note that a normalization factor is applied automatically. In this example, the center-of-mass mapping scheme is used, so the weighting factors for coordinates are defined as the atomic masses.

.. math::

   X_{I,CG} = \frac{ \sum_{i}w_{x,i}x_{i,AA} }{\sum_{i} w_i}


* The `f-weight` list defines the weighting factors to map the CG **forces** in a similar way. The only difference is that there's no automatic normalization when mapping CG forces.

.. math::

   F_{I,CG} = \sum_{i}w_{f,i}f_{i,AA}

In all, the definition of CG types gives two equations that define how to use a
subset of atom coordinates or forces to construct the coordinates or forces for
CG sites, while the `index` list selects the specific atoms that contribute.


CG System
^^^^^^^^^

The `system` section is a list of entries, each of which defines a group of CG sites with a repeated pattern. Each entry in the list must have the following four sections:

* `anchor`: defines the starting index in AA configuration to construct the group of CG sites.
* `repeat`: the number of CG sites defined in this group (e.g., the number of times to repeat the mapping pattern).
* `offset`: the index offset (`anchor`) for each CG site (often the number of
  atoms mapping to the CG bead).
* `sites` or `groups`: the pattern to define the CG sites in this group. In the example of the water box, the `site` section is used. For a more complicated CG topology, the `group` section can be used, which will be described later. Each entry in `sites` is a list of length 2, with the first entry being the site type name and the second entry being the offset in the group where the `site-type` rules are applied.

In the case of the water box, the CG mapping is defined to use three atoms in a
sequence to build a CG site. Therefore, `repeat` is defined as **256** here
for 256 CG sites, while `offset` is defined here as **3**. Please note that the
field `anchor` is define as **0**, so the effective index of anchor atoms during
the mapping are 0, 3, 6, etc.

To map each CG site, the atom index 0 defined by the CG-site type will be
aligned with the anchor atoms. In this example, the CG mapping of WAT will be
aligned with the first atom (index 0) in the AA configuration to build up the
first WAT site. Then, it will slide for 3 atoms, and the atom of index 3 will be
used as the anchor (0 in the CG-site type index) to build up the second WAT
site, and so on. 

Note that there could be more than one scheme of defining a CG mapping.  For
example,

.. code-block::
   :caption: Table 3: Another CG mapping of the water box
   
   site-types:
     WAT:
       index:    [  -1,   0,   1]
       x-weight: [16.0, 1.0, 1.0]
       f-weight: [ 1.0, 1.0, 1.0]
   system:
     - anchor: 1
       repeat: 256
       offset: 3
       sites:
         - [WAT, 0]

This mapping rule is equivalent to the content in Table 1, while the only difference is that the anchor atoms are the indices 1, 4, 7 ... in the AA configurations. However, in the CG-site type, the indices are defined as -1, 0 and 1, so it is still the same groups of atoms, e.g. (0, 1, 2) and (3, 4, 5) are used to construct the CG sites.

Another usage of anchor atoms is for the periodic-boundary-condition (PBC) corrections. When mapping CG coordinates it is assumed that the distances between all atoms that are used for the same CG site and the anchor atom will be corrected to a value less than half of the PBC box. Therefore, the anchors also serve as the anchors for PBC corrections. However, after mapping to CG, the coordinates of CG sites will be corrected using the PBC again to ensure all sites are in the range of the simulation box. Therefore, using different but equivalent CG mapping rules will always give the identical CG trajectory eventually.


Sites Section
^^^^^^^^^^^^^

Each `sites` has a list of sites in the form **[name, offset]**, where the offset updates the anchor atom index when mapping to each site in the list. The example below defines mapping AA methanol into a two-site CG model.

.. code-block::
   :caption: Table 4: CG mapping for two-site methanol
   
   site-types:
     CH3:
       index:    [   0,   1,   2,   3]
       x-weight: [12.0, 1.0, 1.0, 1.0]
       f-weight: [ 1.0, 1.0, 1.0, 1.0]
     OH:
       index:    [   0,   1]
       x-weight: [16.0, 1.0]
       f-weight: [ 1.0, 1.0]
   system:
     - anchor: 0
       repeat: 10
       offset: 6
       sites:
         - [CH3, 0]
         - [OH,  4]

In the mapping group of this example, every six atoms will be mapped to two CG sites, CH3 and OH. Therefore, in the sites section, two CG sites are defined for every six atoms in a group. The first CG site is CH3 and the mapping will be anchored to the atom of index 0, while the mapping of the second CG site OH will be anchored to the atom of index 4.

Multiple groups of mapping can be defined in the system as well. For example,

.. code-block::
   :caption: Table 5: CG mapping for two-site methanol
   
   site-types:
     CH3:
       index:    [   0,   1,   2,   3]
       x-weight: [12.0, 1.0, 1.0, 1.0]
       f-weight: [ 1.0, 1.0, 1.0, 1.0]
     OH:
       index:    [   0,   1]
       x-weight: [16.0, 1.0]
       f-weight: [ 1.0, 1.0]
   system:
     - anchor: 0
       repeat: 10
       offset: 6
       sites:
         - [CH3, 0]
     - anchor: 4
       repeat: 10
       offset: 6
       sites:
         - [OH, 0]

This mapping is equivalent as the previous one, in which all 10 CH3 sites will be mapped first and then all 10 OH sites will be mapped. From the mapping in Table 4, the resulting CG configuration looks like:

.. code-block::
   :caption: CG configuration
   
   ID  Type  X   Y   Z
   1   CH3   ... ... ...
   2   OH    ... ... ...
   3   CH3   ... ... ...
   4   OH    ... ... ...
   5   CH3   ... ... ...
   6   OH    ... ... ...
   ... ..    ... ... ...

While from the mapping in Table 5 the resulting CG configuration looks like:

.. code-block::
   :caption: CG configuration
   
   ID  Type  X   Y   Z
   1   CH3   ... ... ...
   2   CH3   ... ... ...
   3   CH3   ... ... ...
   ... ..    ... ... ...
   11  OH    ... ... ...
   12  OH    ... ... ...
   13  OH    ... ... ...
   
Groups Section
^^^^^^^^^^^^^^

In a mapping group, the `sites` section can be replaced by `groups`, which defines a list of sub mapping groups to be repeated. Therefore, the definition of a mapping group can be recursive. It can largely simplify the mapping definitions for complex models with lots of repeated atom groups. The top-level `system` section can be considered the top-level mapping group. The example below show the CG mapping for an ion-liquid system:

.. admonition:: output

    .. image:: ../_static/cgmap_il.png
        :align: center

The CG mapping can be defined as:

.. code-block::
   :caption: Flattened CG mapping
   
   site-types:
     CH3:
       index:    [   0,   1,   2,   3]
       x-weight: [12.0, 1.0, 1.0, 1.0]
       f-weight: [ 1.0, 1.0, 1.0, 1.0]
     CH2:
       index:    [   0,   1,   2]
       x-weight: [12.0, 1.0, 1.0]
       f-weight: [ 1.0, 1.0, 1.0]
     IMI:
       index:    [   0,    1,   2,    4,    4,   5,    6,   7]
       x-weight: [14.0, 12.0, 1.0, 14.0, 12.0, 1.0, 12.0, 1.0]
       f-weight: [ 1.0,  1.0, 1.0,  1.0,  1.0, 1.0,  1.0, 1.0]
     NO3:
       index:    [   0,    1,    2,    3]
       x-weight: [14.0, 16.0, 16.0, 16.0]
       f-weight: [ 1.0,  1.0,  1.0,  1.0]
   system:
     - anchor: 0
       repeat: 16
       offset: 25
       sites:
         - [CH3,  0]
         - [IMI,  4]
         - [CH2, 12]
         - [CH2, 15]
         - [CH2, 18]
         - [CH3, 21]
     - anchor: 400
       repeat: 16
       offset: 4
       sites:
         - [NO3,  0]

In this example, there are 16 pairs of comprised of a 1-Butyl-3-methylimidazolium cation (25 atoms) and a nitrate anion (4 atoms). The atoms of the cation are mapped to 6 CG sites and atoms in the anion are mapped to 1 CG site. As there are three repeated CH3 sites in the cation CG molecule, the following definition for the cation with recursive mapping groups is equivalent to the previous mapping:

.. code-block::
   :caption: Recursive CG mapping
   
   site-types:
     ...
     ...
   system:
     - anchor: 0
       repeat: 16
       offset: 25
       groups:
         - anchor: 0
           repeat: 1
           offset: 12
           sites:
             - [CH3,  0]
             - [IMI,  4]
         - anchor: 12
           repeat: 3
           offset: 3
           sites:
             - [CH2, 0]
         - anchor: 21
           repeat: 1
           offset: 4
           sites:
             - [CH3,  0]
     - anchor: 400
       repeat: 16
       offset: 4
       sites:
         - [NO3,  0]

Examples
--------

::
    
    cgmap --traj tests/data/methanol_1728_aa.trr \
          --map tests/data/methanol_1728_2s_map.yaml \
          --out cg.trr
