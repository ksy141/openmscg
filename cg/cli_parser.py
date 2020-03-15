
import argparse

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