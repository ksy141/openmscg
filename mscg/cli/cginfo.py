
from pkg_resources import get_distribution
import importlib, inspect, os



def terminal_size():
    import fcntl, termios, struct
    th, tw, hp, wp = struct.unpack('HHHH',
        fcntl.ioctl(0, termios.TIOCGWINSZ,
        struct.pack('HHHH', 0, 0, 0, 0)))
    return tw, th


def main():
    
    # print title
    
    title = "OpenCG Python Package"
    sw = terminal_size()[0]
        
    info = "\n" + " " * ((sw - len(title))//2) + title + "\n" 
    info += " " * ((sw - len(title))//2-1) + "=" * (len(title)+2) + "\n"
    info += "\n\n> Package Information\n---------------------\n\n"
    
    # format metadata
    
    pkg = get_distribution('opencg')
    
    if pkg.has_metadata('PKG-INFO'):
        metadata = list(pkg.get_metadata_lines('PKG-INFO'))
    elif pkg.has_metadata('METADATA'):
        metadata = list(pkg.get_metadata_lines('METADATA'))
    else:
        metadata = []
    
    md = []
    
    import cg
    metadata.append("Location: " + os.path.dirname(cg.__file__))
    
    for row in metadata:
        w = row.split(": ")
        if len(w) == 0:
            continue
        elif len(w) == 1:
            w = ['Description', w[0]]
        else:
            w = [w[0], (": ").join(w[1:])]
        
        if w[0] == 'Requires-Python':
            w[0] = 'Python'
        
        if w[0] != 'Metadata-Version':
            md.append(w)
        
    # print metadata
    
    hw = max([len(row[0]) for row in md])
    
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
    
    # print classes
    
    info += "\n\n> Classes in Package\n--------------------\n\n"
    
    module = importlib.import_module('cg')
    for name, obj in inspect.getmembers(module, inspect.isclass):
        info += " %18s -> %s\n" % (name, str(obj).split("'")[1])
    
    info += "\n\nCongratulations! The package is successfully installed.\n"
    print(info)



if __name__ == '__main__':
    main()




