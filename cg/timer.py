from time import time

class TIMER:
    start = 0.0
    last = 0.0
    recs = {}
    
    @staticmethod
    def reset():
        TIMER.start = time()
        TIMER.last = time()
    
    @staticmethod
    def click(name, f=None):
        now = time()
        
        if name is not None:
            if name not in TIMER.recs:
                TIMER.recs[name] = [0, 0.0]
            
            TIMER.recs[name][0] += 1
            TIMER.recs[name][1] += now - TIMER.last
        
        TIMER.last = now
    
    @staticmethod
    def report(screen = True):
        now = time()
        total = now - TIMER.start
        timed = 0.0
        
        w = 42
        z = []
        
        z.append("=" * w)
        z.append("            Timing Statistics")
        z.append("-" * w)
        z.append(" %-10s %7s %10s %10s" % ("module", "count", "elapsed", "fraction"))
        z.append("-" * w)
        
        for name in TIMER.recs:
            rec = TIMER.recs[name]
            timed += rec[1]
            z.append(" %-10s %7d %10.2f %9.1f%%" % (name[:10], rec[0], rec[1], rec[1]/total * 100))
        
        z.append("-" * w)
        other = total - timed
        z.append(" %-10s %7s %10.2f %9.1f%%" % ("other", "", other, other/total * 100))       
        z.append(" %-10s %7s %10.2f %9.1f%%" % ("total", "", total, 100.0))
        z.append("=" * w)
        
        if screen:
            print("\n".join(z))
        
        return z
        