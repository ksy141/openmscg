import pytest
from cg.cli import cgfm


def test_lj_fluid(datafile):
    
    coeffs = cgfm.main(
        top    = datafile("unary_lj_fluid.top") + ",cgtop",
        traj   = datafile("unary_lj_fluid.lammpstrj,0,1,20"),
        cut    = 2.50,
        pair   = ['1,1,0.9,0.1'],
        verbos = 0,
        save   = 'return'
    )
    
    benchmark = [138.9631694877613, 91.82548595198736, 37.170941251396, 8.274239510050489, -1.2044398613034832, -2.7183788478145345, -2.2683098367441676, -1.5965713385228213, -1.0818547831763554, -0.717259467246318, -0.47796955504588584, -0.3158498732183334, -0.20987270266673158, -0.1361859520093368, -0.08776448922137801, -0.052583077806601486, -0.029078785840772624, -0.014884869849257272, -0.0069221125135834, -0.0021800358760337927, -6.480208702052386e-05]
    
    print("")
    
    for i in range(15):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01



def test_methanol(datafile):
    
    coeffs = cgfm.main(
        top    = datafile("methanol_1728.data") + ",lammps",
        traj   = datafile("methanol_1728.trr"),
        names  = 'MeOH',
        cut    = 8.0,
        pair   = ['MeOH,MeOH,2.8,0.2'],
        verbos = 0,
        save   = 'return'
    )
    
    benchmark = [26.140070724125025, 22.835516255970717, 13.998867238944893, 7.986036757782764, 2.134282485108854, 0.9312967177094646, -1.6252479124254258, 0.015689738206772927, 2.4281807869432646, 0.7596633873710629, 0.7960661642908495, 0.22000532516396565, 0.20661079190148354, -0.05384584359419333, -0.01852867101242014, -0.05076575616135334, -0.023898311446443943, -0.07624874285437323, -0.10042332119456278, -0.047322846328201815, -0.015000554061308858, 0.11587749452651008, 0.05578486627931284, 0.08851018832521618, 0.05291285893352814, 0.0495059268163878, 0.03666081489725824, 0.026932966939835788, 0.0034863381287041344, 0.009183442961512139, 0.015820992008369]
    
    print("")
    
    for i in range(20):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i, benchmark[i], coeffs[i], diff))
        assert abs(diff)<0.01
