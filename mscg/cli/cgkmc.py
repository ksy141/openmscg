import os, numpy as np, matplotlib.pyplot as plt

global colorlist, vmd_color_ids
colorlist = ["royalblue", "crimson", "mediumorchid", "lightseagreen",
"greenyellow", "slateblue", "olive", "yellow", "orange", "limegreen", "lightsteelblue", "pink", "darkblue", 
"deepskyblue", "darkkhaki", "hotpink", "yellowgreen", "maroon", "rosybrown", "teal"]
vmd_color_ids = [0, 1, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15]

def similarity(list_0, list_1):
    similar_elements = []
    for item in list_1:
        if item in list_0:
            similar_elements.append(item)
    return len(similar_elements)**2 / len(list_0) / len(list_1)

class CG_cluster(object):
    def __init__(self) -> None:
        self.all_atom_num = 0
        self.frame_num = 0
        self.frames = []
        self.alpha = 0.1 # Sets an acceptance rate for chi_square elevating steps.
        self.beta = 1 # Controls the contribution of spatial closeness to the overall chi_square.
        self.gamma = 0 # Controls the contribution of sequence intactness to the overall chi_square.
        self.mode = 0 # 0 for stocastic and 1 for sequence cutter.

    def read_parm_from_psf(self, filename):
        self.parmfile = filename
        with open(filename) as f:
            print("Reading parameters from psf file", filename)
            while True:
                a = f.readline()
                if "NATOM" in a:
                    break
            
            self.atom_labels = []
            self.residue_labels = []
            self.monomer_labels = []
            self.residue_ids = []
            self.charges = []
            self.masses = []
            self.CA_atom_ids = []

            while True:
                a = f.readline()
                content = a.strip().split()
                if len(content) != 9: break
                self.atom_labels.append(content[4])
                if content[4] == "CA":
                    self.CA_atom_ids.append(int(content[0])-1)
                self.residue_labels.append(content[3])
                self.monomer_labels.append(content[1])
                self.residue_ids.append(int(content[2]))
                self.charges.append(float(content[-3]))
                self.masses.append(float(content[-2]))

    def read_frames_from_rst7(self, filename, only_one_residue = False, only_read_CA = True, max_frame_num = 0):
        self.trajfile = filename
        print("Reading frames from rst7 file", filename)
        if only_one_residue == True:
            self.stop_point = 0
            while self.monomer_labels[self.stop_point+1] == self.monomer_labels[self.stop_point]:
                self.stop_point += 1
        
        with open(filename) as f:
            line_num = len(f.readlines())

        with open(filename) as f:        
            frame = []
            line_id = 0
            
            while line_id < line_num:
                line = f.readline()
                content = line.strip().split()
                if line_id == 0:
                    pass
                elif line_id == 1:
                    self.all_atom_num = int(content[0])
                    if only_one_residue == False: 
                        self.stop_point = self.all_atom_num
                    lines_per_frame = np.ceil(self.all_atom_num/2) + 2
                    self.frame_num += 1
                    if self.frame_num % 200 == 0: print("Reading information from frame", self.frame_num)
                elif len(content) == 1:
                    # This marks the start of every frame. We must clear the previous frames.
                    frame = []
                    self.frame_num += 1
                    if max_frame_num > 0 and self.frame_num > max_frame_num: break
                    if self.frame_num % 200 == 0: print("Reading information from frame", self.frame_num)
                elif (self.frame_num-1)*lines_per_frame+1 < line_id < (self.frame_num)*lines_per_frame:
                    current_atom_id = 2*(line_id - (self.frame_num-1)*lines_per_frame+1)
                    if current_atom_id > self.stop_point: pass
                    elif len(content) == 6:
                        frame.append(np.asarray(content[:3], dtype = np.float64))
                        frame.append(np.asarray(content[3:6], dtype = np.float64))
                    elif len(content) == 3:
                        frame.append(np.asarray(content[:3], dtype = np.float64))
                elif line_id == self.frame_num * lines_per_frame:
                    frame = np.asarray(frame)
                    if only_read_CA == True:
                        CA_atom_ids = []
                        for i in self.CA_atom_ids:
                            if i > self.stop_point: break
                            CA_atom_ids.append(i)
                        self.CA_atom_ids_before_stop = np.asarray(CA_atom_ids)
                        self.frames.append(frame[self.CA_atom_ids_before_stop])
                    else:
                        self.frames.append(frame)
                    # Every line contains two coordinates.
                    # The end line is marked with the box size and angles.

                line_id += 1

        self.frames = np.asarray(self.frames)
        self.atom_num = len(self.frames[0])
        self.eq_positions = np.mean(self.frames, axis=0)
        self.fluctuation_traj = self.frames - np.asarray([self.eq_positions for k in range(self.frame_num)])

    def initialize_labels(self, num=12, seed=None):
        """
        Initializes labels for the sequence.
        ---
        num = target label number. Set to 12 for default.
        random = assign initial labels for each residue randomly.

        """

        self.label_num = num
        self.labels = np.zeros(self.atom_num)
        self.reverse_labels = [[] for k in range(num)]
        batch_num = int(np.ceil(self.atom_num / num))
        if self.mode == 0:
            rng = np.random.default_rng(seed=seed)
            self.labels = rng.integers(low=0, high=self.label_num, size=self.atom_num)
            for i in range(len(self.labels)):
                self.reverse_labels[self.labels[i]].append(i)
        elif self.mode == 1:
            rng = np.random.default_rng(seed=seed)
            cutters = rng.integers(low=1, high=(self.atom_num-1), size=(self.label_num-1)*10)
            cutters = np.unique(cutters)
            np.random.shuffle(cutters)
            self.cutters = sorted(list(cutters[:self.label_num-1]))
            self.cutters.append(self.atom_num)

            self.reverse_labels[0] = list(range(0, self.cutters[0]))
            for i in range(len(self.cutters)-1):
                self.labels[self.cutters[i]:self.cutters[i+1]] = i+1
                self.reverse_labels[i+1] = list(range(self.cutters[i], self.cutters[i+1]))

        self.labels = list(np.asarray(self.labels, dtype=int))
        for i in range(num):
            self.reverse_labels[i].sort()

    def read_labels_from_file(self, filename, only_return=False):
        reverse_labels = []
        labels = []
        with open(filename) as f:
            aa = f.readlines()
            for line in aa:
                lst = eval(line)
                if type(lst) == list:
                    reverse_labels.append(lst)
                else:
                    labels.append(lst)
        label_list = np.zeros(self.atom_num)
        for reverse_label in reverse_labels:
            label = reverse_labels.index(reverse_label)
            for atom_id in reverse_label:
                label_list[atom_id] = label
        
        if only_return == True:
            return list(np.asarray(label_list, dtype=int)), reverse_labels
        else:
            self.labels, self.reverse_labels = list(np.asarray(label_list, dtype=int)), reverse_labels
            self.label_num = len(labels)

    def update_label(self, atom_id, old_label, new_label):
        self.labels[atom_id] = new_label
        self.reverse_labels[old_label].remove(atom_id)
        self.reverse_labels[old_label].sort()
        self.reverse_labels[new_label].append(atom_id)
        self.reverse_labels[new_label].sort()

    def calculate_rmsf(self):
        try: os.mkdir('cache_0728')
        except: pass
        print("Calculating RMSF list for all atoms:")
        self.rmsf_list = np.mean(np.linalg.norm(self.fluctuation_traj, axis=2)**2, axis=0) ** (1/2)
        self.msf_list = np.mean(np.linalg.norm(self.fluctuation_traj, axis=2)**2, axis=0)
        np.save("cache_0728/rmsf.npy", self.rmsf_list)
            
    def calculate_chi_square(self):

        self.num_cuts = []
        self.label_centers = []
        self.label_center_fluc = np.zeros((self.frame_num, self.label_num, 3))

        self.chi_list = np.zeros((3, self.label_num))

        for i in range(self.label_num):
            gom = np.mean(self.eq_positions[self.reverse_labels[i]], axis=0, keepdims=True)
            gom_flat = np.mean(self.eq_positions[self.reverse_labels[i]], axis=0)
            self.label_centers.append(gom_flat)
            gom_fluc_flat = np.mean(self.fluctuation_traj[:,self.reverse_labels[i],:], axis=1)
            gom_fluc = np.mean(self.fluctuation_traj[:,self.reverse_labels[i],:], axis=1, keepdims=True)
            self.label_center_fluc[:,i,:] = gom_fluc_flat
            num_cuts = 0
            for j in range(len(self.reverse_labels[i])-1):
                if self.reverse_labels[i][j+1] - self.reverse_labels[i][j] > 1: num_cuts += 1
            self.num_cuts.append(num_cuts)
            chi_1 = np.mean(np.sum(np.linalg.norm(self.fluctuation_traj[:,self.reverse_labels[i],:] - gom_fluc, ord=2, axis=2)**2, axis=1))
            chi_2 = np.sum(np.linalg.norm(self.eq_positions[self.reverse_labels[i]] - gom, ord=2, axis=1)**2)
            chi_3 = num_cuts ** 2
            self.chi_list[0][i] = chi_1
            self.chi_list[1][i] = chi_2
            self.chi_list[2][i] = chi_3

        [self.chi_1, self.chi_2, self.chi_3] = np.sum(self.chi_list, axis=1)

        chi_square = self.chi_1 + (self.beta * self.chi_2) / self.atom_num + self.gamma * self.chi_3
        chi_square = chi_square / 3  
        self.chi_square = chi_square
        self.label_centers = np.asarray(self.label_centers)

    def update_chi_square(self, atom_id, old_label, new_label):

        labels = self.labels.copy()
        new_num_cuts = self.num_cuts.copy()
        chi_list = self.chi_list.copy()
        labels[atom_id] = new_label

        for i in (old_label, new_label):
            indices = self.reverse_labels[i].copy()
            if i == old_label: indices.remove(atom_id)
            else: indices.append(atom_id)
            indices.sort()

            gom = np.mean(self.eq_positions[indices], axis=0, keepdims=True)
            gom_fluc = np.mean(self.fluctuation_traj[:,indices,:], axis=1, keepdims=True)
            num_cuts = 0
            for j in range(len(indices)-1):
                if indices[j+1] - indices[j] > 1: num_cuts += 1
            new_num_cuts[i] = num_cuts
            chi_1 = np.mean(np.sum(np.linalg.norm(self.fluctuation_traj[:,indices,:] - gom_fluc, ord=2, axis=2)**2, axis=1))
            chi_2 = np.sum(np.linalg.norm(self.eq_positions[indices] - gom, ord=2, axis=1)**2)
            chi_3 = num_cuts ** 2
            chi_list[0][i] = chi_1
            chi_list[1][i] = chi_2
            chi_list[2][i] = chi_3

        [c1, c2, c3] = np.sum(chi_list, axis=1)

        chi_square = c1 + (self.beta * c2) / self.atom_num + self.gamma * c3
        return chi_square / 3

    def optimize_chi_square(self, cluster_converge_steps = 100, max_steps=10000, output_interval = 100, filename = None):

        if filename == None: file_str = "%.0f" % (time() % (5*86400))
        else: file_str = str(filename)
        if output_interval != 0: chi_square_traj_file = open("chi_square_traj_%s.txt" % file_str, "w")
        self.chi_square_traj = []

        for step in range(cluster_converge_steps):
            self.calculate_chi_square()
            self.chi_square_traj.append(self.chi_square)

            if output_interval != 0 and step % output_interval == 0:
                print("Step = ", step, "\tChi_square = %.6f = %.6f + %.6f + %.6f" % (self.chi_square, self.chi_1, self.chi_2, self.chi_3))
                print("%d\t%.6f" % (step, self.chi_square), file=chi_square_traj_file)

            if len(self.chi_square_traj) > 21 and np.var(self.chi_square_traj[-20:]) < 1e-6:
                print("Initial optimization complete. Steps = %d, Chi_square = %.4f" % (step+1, self.chi_square))
                break

            for i in range(self.atom_num):
                old_label = self.labels[i].copy()
                term_1 = self.beta * np.linalg.norm(self.label_centers - self.eq_positions[i], ord=2, axis=1)
                term_2 = np.mean(np.linalg.norm(self.label_center_fluc - self.fluctuation_traj[:,i:(i+1),:], axis=2)**2, axis=0)
                dist_list = term_1 + term_2
                new_label = np.argmin(dist_list)
                if len(self.reverse_labels[old_label]) > 1: 
                    self.update_label(i, old_label, new_label)

        step_count = step

        rng = np.random.default_rng()
        successes = []
        trial_CA_atoms = rng.integers(low=0, high=self.atom_num, size=max_steps)
        domain_seeds = np.asarray([(np.asarray(self.labels)+k)%self.label_num for k in range(max_steps//self.atom_num+1)], dtype=int).flatten()
        #domain_seeds = rng.integers(low=0, high=self.label_num, size=max_steps)

        for step in range(max_steps):
            reverse_labels = self.reverse_labels.copy()
            if output_interval != 0 and (step) % output_interval == 0:
                print("Step = ", step, "\tChi_square = %.6f = %.6f + %.6f + %.6f" % (self.chi_square, self.chi_1, self.chi_2, self.chi_3))
                print("%d\t%.6f" % (step, self.chi_square), file=chi_square_traj_file)

            if len(successes) > 5001 and sum(successes[-5000:]) < 1e-6:
                print("Second phase optimization complete. Steps = %d, Chi_square = %.4f" % (step+1, self.chi_square))
                break

            #trial_atom = trial_CA_atoms[step]
            trial_atom = step%self.atom_num
            old_label = self.labels[trial_atom]
            if len(self.reverse_labels[old_label]) <= 1: pass
            
            new_label = domain_seeds[step].copy()
            if new_label == old_label:
                new_label = (new_label+1) % self.label_num

            new_chi_square = self.update_chi_square(trial_atom, old_label, new_label)
            if new_chi_square < self.chi_square:
                self.update_label(trial_atom, old_label, new_label)
                self.calculate_chi_square()
                successes.append(1)
            else: 
                successes.append(0)

            self.chi_square_traj.append(self.chi_square)

        self.write_labels(filename)
        self.write_labels_clean(filename)

        if output_interval != 0: chi_square_traj_file.close()

    def align_labels(self, standard_reverse_labels, this_reverse_labels, only_return=False):
        assert len(standard_reverse_labels) == len(this_reverse_labels)
        new_reverse_labels = [0 for k in range(len(this_reverse_labels))]
        new_indices = []
        lst = list(range(len(standard_reverse_labels)))
        for k in range(len(this_reverse_labels)):
            similarity_list = -np.ones(len(this_reverse_labels))
            for j in lst:
                similarity_list[j] = similarity(standard_reverse_labels[j], this_reverse_labels[k])
            new_index = np.argmax(similarity_list)
            new_indices.append(new_index)
            lst.remove(new_index)
            new_reverse_labels[new_index] = this_reverse_labels[k]
        if only_return == True:
            return new_indices, new_reverse_labels
        else:
            self.reverse_labels = new_reverse_labels
            for k in range(self.atom_num):
                old_label = self.labels[k].copy()
                self.labels[k] = new_indices[old_label]

    def align_labels_from_file(self, standard_set_filename, this_set_filename, only_return=False):
        standard_labels, standard_reverse_labels = self.read_labels_from_file(standard_set_filename, only_return=True)
        this_labels, this_reverse_labels = self.read_labels_from_file(this_set_filename, only_return=True)
        new_indices, new_reverse_labels = self.align_labels(standard_reverse_labels, this_reverse_labels)
        if only_return == True:
            return new_indices, new_reverse_labels
        else:
            self.reverse_labels = new_reverse_labels
            for k in range(self.atom_num):
                old_label = self.labels[k].copy()
                self.labels[k] = new_indices[old_label]

    def write_labels(self, filename=None):
        if filename == None: file_str = "%.0f" % (time() % (5*86400))
        else:file_str = str(filename)
        with open("result_"+file_str+".txt", 'w') as f_clean:
            for k in range(self.label_num):
                print(k, file=f_clean)
                print(self.reverse_labels[k], file=f_clean)

    def write_labels_clean(self, filename=None):
        if filename == None: file_str = "%.0f" % (time() % (5*86400))
        else: file_str = str(filename)
        with open("result_clean_"+file_str+".txt", 'w') as f_clean:
            for k in range(self.label_num):
                print(k, file=f_clean)
                out_str = "[["+str(self.reverse_labels[k][0])+","
                for i in range(len(self.reverse_labels[k])-1):
                    if self.reverse_labels[k][i+1] - self.reverse_labels[k][i] > 1.01:
                        out_str += str(self.reverse_labels[k][i])+"],["+str(self.reverse_labels[k][i+1])+","
                out_str += str(self.reverse_labels[k][-1]) + "]]"
                print(out_str, file=f_clean)

    def visualize_labels(self, parmfile:str = None, trajfile:str = None, starter:int=0):
        if parmfile == None: parmfile = self.parmfile
        if trajfile == None: trajfile = self.trajfile

        with open("vmd_script.in", "w") as vmd:
            print("""mol new {0} 
mol addfile {1}
mol color name
mol representation NewCartoon
mol material Opaque""".format(parmfile, trajfile), file=vmd)
            for k in range(len(self.reverse_labels)):
                print("mol addrep", 0, file=vmd)
                reverse_labels = list(starter + np.asarray(self.reverse_labels[k], dtype=int))
                select_range = "residue %s" % " ".join(list(map(str, reverse_labels)))
                print("mol modselect", k+1, 0, select_range, file=vmd)
                print("mol modcolor", k+1, 0, "ColorID", vmd_color_ids[k], file=vmd)
            print("""mol delrep 0 0
color Display Background white
color change rgb 8 white 1.000000 1.000000 1.000000""", file=vmd)
        os.system("vmd -e vmd_script.in")

    def plot_labels(self, reverse_labels_1 = None, reverse_labels_2 = None, ref_reverse_labels = None, filename = None):
        if filename == None: file_str = "%.0f" % (time() % (5*86400))
        else: file_str = str(filename)

        fig, ax = plt.subplots(dpi = 200, figsize = (12, 4), tight_layout=True)
        ax.plot(range(len(self.rmsf_list)), self.rmsf_list, color = 'black', alpha=0.2)

        if reverse_labels_1 == None: 
            reverse_labels_1 = self.reverse_labels

        for lst in reverse_labels_1:
            k = reverse_labels_1.index(lst)
            ax.scatter(lst, self.rmsf_list[lst], color = colorlist[k], s = 9, label = "CG site #%d" % (k+1))

        if reverse_labels_2 != None:
            for lst in reverse_labels_2:
                k = reverse_labels_2.index(lst)
                ax.scatter(lst, self.rmsf_list[lst]+0.15, color = colorlist[k], s = 9)
                
        if ref_reverse_labels != None:
            for k in range(len(ref_reverse_labels)):
                out_str = "[["+str(ref_reverse_labels[k][0])+","
                for i in range(len(ref_reverse_labels[k])-1):
                    if ref_reverse_labels[k][i+1] - ref_reverse_labels[k][i] > 1.01:
                        out_str += str(ref_reverse_labels[k][i])+"],["+str(ref_reverse_labels[k][i+1])+","
                out_str += str(ref_reverse_labels[k][-1]) + "]]"
                label = eval(out_str)
                for sub_label in label:
                    ax.fill_between(range(sub_label[0], sub_label[1]+1), 0, self.rmsf_list[sub_label[0]:sub_label[1]+1], alpha=0.3, color = colorlist[k])
        else:
            for k in range(len(reverse_labels_1)):
                out_str = "[["+str(reverse_labels_1[k][0])+","
                for i in range(len(reverse_labels_1[k])-1):
                    if reverse_labels_1[k][i+1] - reverse_labels_1[k][i] > 1.01:
                        out_str += str(reverse_labels_1[k][i])+"],["+str(reverse_labels_1[k][i+1])+","
                out_str += str(reverse_labels_1[k][-1]) + "]]"
                label = eval(out_str)
                for sub_label in label:
                    ax.fill_between(range(sub_label[0], sub_label[1]+1), 0, self.rmsf_list[sub_label[0]:sub_label[1]+1], alpha=0.3, color = colorlist[k])
                
        ax.set_xlim(-10, len(self.rmsf_list)+10)
        ax.set_ylim(0, max(self.rmsf_list)+1)
        ax.set_xlabel('CA Atom Index')
        ax.set_ylabel('RMSF (angstrom)')
        ax.set_title("Cluster result")
        if len(reverse_labels_1) < 12: ax.legend()
        fig.savefig("result_plot_%s.png" % file_str, dpi=400)

if __name__ == '__main__':
    main()
