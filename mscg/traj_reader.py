from .trajectory import Trajectory
from .verbose import screen
from time import time


class TrajReader:
    
    def __init__(self, file, skip, every, frames):
        
        segs = file.split(".")
        suffix = segs[-1] if len(segs)>1 else ""
        
        self.traj   = None
        self.file   = file
        self.suffix = suffix
        self.skip   = skip
        self.every  = every
        self.frames = frames
        self.nread  = 0
        
    def close(self):
        del(self.traj)
        self.traj = None
            
    def next_frame(self):
        
        if self.traj is None:
            self.nread = 0
            self.traj = Trajectory(self.file, self.suffix)
            for i in range(self.skip):
                self.traj.read_frame()
        
        if self.frames>0 and self.nread>=self.frames:
            self.close
            return False
        
        for i in range(self.every-1):
            self.traj.read_frame()
        
        if self.traj.read_frame():
            self.nread += 1
            return True
        else:
            self.close()
            return False



class TrajBatch:
    def __init__(self, readers, natoms = None, cut = None):
        self.readers = readers
        self.natoms = natoms
        self.cut = cut
        
    def __iter__(self):
        self.reader_iter = iter(self.readers)
        self.reader = None
        return self
    
    def __next__(self):
        if self.reader is None:
            success = False
        else:
            success = self.reader.next_frame()
            
            if not success:
                elapsed = time() - self.timer_start
                msg = " -> Processed %d frames. Elapsed: %0.0f secs." % (self.reader.nread, elapsed)
                screen.info(('\r%s' % (msg)) + " " * 30)
        
        while not success:
            self.reader = next(self.reader_iter, None)
            
            if self.reader is None:
                break
            else:
                success = self.reader.next_frame()
                
                if success:
                    if (self.natoms is not None) and (self.natoms != self.reader.traj.natoms):
                        screen.fatal("Inconsistent number of atoms between topology (%d) and trajectory (%d)." % (self.natoms, self.reader.traj.natoms))

                    if self.cut is not None:
                        cut2 = self.cut * 2.0
                        box = self.reader.traj.box
                        
                        if box[0]<cut2 or box[1]<cut2 or box[2]<cut2:
                            screen.fatal("Incorrect cut-off for the trajectory: cut-off (%f) must be smaller than half of the box dimentions (%s)" % (self.cut, str(box)))
                                        
                    screen.info("Process trajectory: " + self.reader.file)
                    self.timer_start = self.timer_last = time()

        if success:
            now = time()
            
            if now - self.timer_last > 1.0:
                self.timer_last = now
                elapsed = now - self.timer_start

                if self.reader.frames>0:
                    remained = (now - start) / self.reader.nread * (self.reader.frames - self.reader.nread)
                    msg = " -> Processed %d of %d frames. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (self.reader.nread, self.reader.frames, elapsed, remained)
                else:
                    msg = " -> Processed %d frames. Elapsed: %0.0f secs ..." % (self.reader.nread, elapsed)

                screen.info('\r%s' % (msg), end="")

            return self.reader
        else:
            raise StopIteration
        
            
        
        
        
        
        
        
        
        
        
    
