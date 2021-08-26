
import argparse



def parse_kwargs_line(line):
    kwargs = {}
    
    for k, v in [seg.split('=') for seg in line.split(',')]:
        if v.replace('.','').isnumeric():
            kwargs[k] = float(v)
        elif v.lower() == 'true':
            kwargs[k] = True
        elif v.lower() == 'false':
            kwargs[k] = False
        else:
            kwargs[k] = v
    
    return kwargs

class CLIParser(argparse.ArgumentParser):
    
    def convert_arg_line_to_args(self, arg_line):
        if arg_line.strip()[:1] == '#':
            return []
        
        return arg_line.strip().split()
    
    def parse_inline_args(self, *args, **kwargs):
        
        arg_list = list(args).copy()
        
        for k,v in kwargs.items():
            
            if type(v) == bool:
                if v == True:
                    arg_list.append('--' + k)
            elif type(v) == list:
                for one in v:
                    arg_list.append('--' + k)
                    arg_list.append(str(one))
            else:
                arg_list.append('--' + k)
                arg_list.append(str(v))
        
        return self.parse_args(arg_list)



from .topology import Topology

class TopAction(argparse.Action):
    
    help = "topology file"
    
    def __call__(self, parser, namespace, values, option_string=None):
        
        if type(values) != str:
            raise ValueError("incorrect format of value for option --top")
        
        top = Topology.read_file(values)
        setattr(namespace, self.dest, top)


from .traj_reader import TrajReader

class TrajReaderAction(argparse.Action):
    
    help = "reader for a trajectory file, multiple fields separated by commas, the first field is the file name, while others define the skip, every and frames (default args: file,skip=0,every=1,frames=0)"
    
    def __call__(self, parser, namespace, values, option_string=None):
        
        if type(values) != str:
            raise ValueError("incorrect format of value for option --traj")
        
        skip, every, frames = 0, 1, 0
        segs = values.split(",")
        
        for i in range(1, len(segs)):
            w = [seg.strip() for seg in segs[i].split("=")]
            
            if len(w)==1:
                raise ValueError("incorrect format of value for option --traj: " + segs[i])
            elif w[0] == 'skip':
                skip = int(w[1])
            elif w[0] == 'every':
                every = int(w[1])
            elif w[0] == 'frames':
                frames = int(w[1])
            else:
                raise ValueError("incorrect format of value for option --traj: " + segs[i])
        
        reader = TrajReader(segs[0], skip, every, frames)
        getattr(namespace, self.dest).append(reader)
        
        
        
