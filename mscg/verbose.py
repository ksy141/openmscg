
class screen:
    
    # verbose level: 0 - silent, 1 - info, 2 - debug
    
    verbose = 1
        
    @staticmethod
    def show(msg, prefix="", end="\n"):
        if type(msg) == str:
            msg = [msg]
        
        print("\n".join([prefix + str(s) for s in msg]), end=end)
    
    @staticmethod
    def error(msg):
        screen.show(msg, "Error: ")
        
    @staticmethod
    def fatal(msg):
        screen.error(msg)
        exit()
    
    @staticmethod
    def info(msg, end="\n"):
        if screen.verbose>0:
            screen.show(msg, "", end)
    
    @staticmethod
    def debug(msg, end="\n"):
        if screen.verbose>1:        
            screen.show(msg, "debug: ", end)
