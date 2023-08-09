import numpy as np
f=open('protein_1microsec.pdb','r')
resid = []
element = []
lines = f.readlines()
#print(type(len(lines)))
for i in range(5,len(lines)-2):
    entries = lines[i].split()
    resid.append(entries[5])
    element.append(entries[11])#11
#print(len(resid))
mass = np.zeros(len(element))
for i in range(len(element)):
    if element[i] == 'C':
        mass[i] = 12.0
    if element[i] == 'H':
        mass[i] = 1.0
    if element[i] == 'O':
        mass[i] = 16.0
    if element[i] == 'N':
        mass[i] = 14.0
    if element[i] == 'S':
        mass[i] = 32.0
atoms_per_AA = []
CA_count = 1
num= 1
for i in range(len(resid)-1):
    if resid[i] != resid[i+1]:
        atoms_per_AA.append(num)
        num = 1
        CA_count += 1
    else:
        num += 1
atoms_per_AA.append(num)        
print(atoms_per_AA) #to be used in the sites subsection in system
#print(len(atoms_per_AA))
#print(CA_count)
mapping = np.loadtxt('map.txt')
mapping_beads = len(mapping)#int(CA_count/mapping)#+1 add +1 for resolutions higher than 1 bead per AA
#print(mapping_beads)
resid_range = np.arange(1,CA_count+1)
#print(mapping_beads)
#print(len(resid_range))
f.close()
cg_site_name = []
for site in range(mapping_beads):
    cg_site_name.append("B"+str(site+1))
filename='resolution_edcg.yaml'
f=open(filename,'w')
f.write("site-types: \n")
f.close()
#atoms_per_site = np.zeros(mapping_beads)
#for i in range(0,mapping_beads):
#    print(atoms_per_AA[int(mapping[i,0]):int(mapping[i,1])+1])

atoms_per_site = []
for i in range(0,mapping_beads):
    if (int(mapping[i,0]) - int(mapping[i,1])) != 0:
        atoms_per_site.append(sum(atoms_per_AA[int(mapping[i,0]):int(mapping[i,1])+1]))
    else:
        atoms_per_site.append(atoms_per_AA[int(mapping[i,0])])
print(atoms_per_site)

f=open(filename,'a')
for i in range(mapping_beads):
    first_index = 0
    f.write(f"        {cg_site_name[i]}: \n")
    index_range = atoms_per_site[i]
    array = list(np.arange(index_range))
    array_fweight = list(np.ones(index_range))
    array_mass = list(mass[first_index:atoms_per_site[i]])
    f.write(f"            index: {array} \n")
    f.write(f"            x-weight: {array_mass} \n")
    f.write(f"            f-weight: {array_fweight} \n")
    first_index += atoms_per_site[i]

f.write("system: \n")
f.write("      - anchor: 0 \n")
f.write("        repeat: 1 \n")
f.write(f"        offset: {len(resid)} \n")
f.write(f"        sites: \n")
count = 0
for i in range(len(cg_site_name)):
    if i == 0:
        f.write(f"            - [{cg_site_name[i]},0] \n")
    else:
        count += atoms_per_site[i-1]
#        #print(count)
        f.write(f"            - [{cg_site_name[i]},{count}] \n")
f.close()

