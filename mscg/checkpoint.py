#

import os, sys, socket, getpass, datetime
import pickle
import __main__

class Checkpoint:
    
    def __init__(self, filename, script):
        
        self.data = {
            'script': script,
            'command': " ".join(sys.argv),
            'os': sys.platform,
            'host': socket.gethostname(),
            'user': getpass.getuser(),
            'path': os.getcwd(),
            'time': str(datetime.datetime.now()),
            'file': filename + '.p',
        }
    
    def update(self, data):
        self.data.update(data)
        return self
    
    def dump(self):
        pickle.dump(self.data, open(self.data['file'], 'wb'))
