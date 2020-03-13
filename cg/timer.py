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
    def report():
        now = time()
        total = now - TIMER.start
        timed = 0.0
        
        w = 42
        print("=" * w)
        print("            Timing Statistics")
        print("-" * w)
        print(" %-10s %7s %10s %10s" % ("module", "count", "elapsed", "fraction"))
        print("-" * w)
        
        for name in TIMER.recs:
            rec = TIMER.recs[name]
            timed += rec[1]
            print(" %-10s %7d %10.2f %9.1f%%" % (name[:10], rec[0], rec[1], rec[1]/total * 100))
        
        print("-" * w)
        other = total - timed
        print(" %-10s %7s %10.2f %9.1f%%" % ("other", "", other, other/total * 100))       
        print(" %-10s %7s %10.2f %9.1f%%" % ("total", "", total, 100.0))
        print("=" * w)