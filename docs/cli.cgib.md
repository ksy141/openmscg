# Command: cgib

The script to calculate distributions of pairs, bonds, angles ... It can also calculate the potential functions from the inversed-boltzmann method.

---

### Arguments:


--skip

--every

--frames

--top

--names

--cut

--temp

--pair

--bond

--angle 


### Example:

ibcg
--skip 0
--every 10
--frames 100
--top lammps
--names CH3,OH
--cut 10.0
--temp 298.15
--pair CH3,CH3,2.8,10.0,200
--pair CH3,OH,2.8,10.0,200
--pair OH,OH,2.5,10.0,200
--bond CH3,OH,1.35,1.65,60
