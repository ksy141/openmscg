
from pkg_resources import get_distribution
import importlib, inspect



def terminal_size():
    import fcntl, termios, struct
    th, tw, hp, wp = struct.unpack('HHHH',
        fcntl.ioctl(0, termios.TIOCGWINSZ,
        struct.pack('HHHH', 0, 0, 0, 0)))
    return tw, th


def main():
    
    title = "OpenCG Python Package"
    
    pkg = get_distribution('opencg')
    md = [row.split(": ") for row in pkg.get_metadata('PKG-INFO').strip().split("\n")[1:]]
    hw = max([len(row[0]) for row in md])
    sw = terminal_size()[0]
    
    info = "\n" + " " * ((sw - len(title))//2) + title + "\n" 
    info += " " * ((sw - len(title))//2-1) + "=" * (len(title)+2) + "\n"
    info += "\n\n> Package Information\n---------------------\n\n"
    
    for row in md:
        h = ("%" + str(hw) +"s: ") % (str(row[0]))
        info += h
        w = 0
        
        for word in row[1].split(" "):
            word += " "
            
            if w + len(word) > sw - hw - 2:
                info += "\n" + " " * (hw + 2)
                w = 0
            
            info += word
            w += len(word)
            
        info += "\n"
    
    info += "\n\n> Classes in Package\n--------------------\n\n"
    
    module = importlib.import_module('cg')
    for name, obj in inspect.getmembers(module, inspect.isclass):
        info += " %18s -> %s\n" % (name, str(obj).split("'")[1])
    
    info += "\n\nCongratulations! The package is successfully installed.\n"
    print(info)



if __name__ == '__main__':
    main()




