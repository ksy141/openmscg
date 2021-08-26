"""
Examples with Reader Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Besides directly calling the low-level APIs in the `Trajectory` class,
the package also provides a high-level `TrajReader` class to simplify 
the reading process of a group of trajectory files with flexible 
controling parameters::

    >>> from mscg import TrajReader
    >>> reader = TrajReader('tests/data/methanol_1728_cg.trr', skip=50, every=10, frames=25)
    >>> nread = 0
    >>> while reader.next_frame():
    ...     nread += 1
    ...
    >>> nread
    25

To read a group of trajectories consecutively, the package also provide
an `iterator` class for the batch processing::
    
    >>> from mscg import TrajBatch
    >>> for reader in TrajBatch(readers):
    ...     if reader.nread == 1:
    ...         pass # do something here for the first frame in a file
    ...
    ...     pass # other operations for the frames
    ...

For the development of a CLI script, a customized argument parser is provided
to accept multiple trajectory arguments for processing. This is widely used by
the CLI commands included in this package:
    
.. code-block:: python
    
    from mscg import *
    
    parser = CLIParser(description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        fromfile_prefix_chars='@', add_help=False)
    
    parser.add_argument("--traj", metavar='file[,args]',
        action=TrajReaderAction, help=TrajReaderAction.help, default=[])
    ...
    
    args = parser.parse_args()
    ...
    
    for reader in TrajBatch(args.traj):
        pass # processing frames in the trajectories specified by --traj
        ...
    ...

More instructions for this argument parser can be found `here <../basics.html#cli-option-for-trajectory>`__

"""

from .trajectory import Trajectory
from .verbose import screen
from time import time


class TrajReader:
    """
    A reader controller class for reading frames from a trajectory file.
    
    Attribute:
        traj : mscg.Trajectory
            a trajectory object controlled by this reader.
    """
    
    def __init__(self, file:str, skip:int=0, every:int=1, frames:int=0):
        """
        Create a reader controller for a trajectory file.
        
        :param file: file name.
        :type file: str
        
        :param skip: number of frames to be skipped at the head, default to 0.
        :type skip: int
        
        :param every: read a frame data for every the number of frames, default to 1.
        :type every: int
        
        :param frames: total number of frames to be read (0 for reading all frames), default to 0.
        :type frames: int
        """
        
        segs = file.split(".")
        suffix = segs[-1] if len(segs)>1 else ""
        
        self.traj   = None
        self.file   = file
        self.suffix = suffix
        self.skip   = skip
        self.every  = every
        self.frames = frames
        self.nread  = 0
    
    def __del__(self):
        self.close()
        
    def close(self):
        """
        Close the reader and release the trajectory object.
        """
        
        if self.traj is not None:
            del(self.traj)
            self.traj = None
            
    def next_frame(self):
        """
        Read the next available frame in the trajectory according
        to the controlling parameters.
        
        :return: return *True* if new frame data is read, else return *False*.
        :rtype: bool
        """
        
        if self.traj is None:
            self.nread = 0
            self.traj = Trajectory(self.file, "r", self.suffix)
            for i in range(self.skip):
                self.traj.read_frame()
        
        if self.frames>0 and self.nread>=self.frames:
            self.close()
            return False
        
        for i in range(self.every-1):
            self.traj.read_frame()
        
        if self.traj.read_frame():
            self.nread += 1
            return True
        else:
            self.close()
            return False
    
    def __iter__(self):
        self.close()
        return self
    
    def __next__(self):
        if self.next_frame():
            return self.traj
        else:
            raise StopIteration



class TrajBatch:
    """
    A batch iterator for multiple trajectory readers. The check enforcements can
    be set that, if the number of atoms or box sizes of the system in a new trajectory
    dosen't match the requirement, an exception will be raised up.
    """
    
    def __init__(self, readers, natoms = None, cut = None):
        """
        Create a batch object.
        
        :param readers: a list of trajectory readers
        :type readers: [mscg.TrajReader]
        
        :param natoms: enforcing certain number atoms in the trajectory if not *None*.
        :type natoms: int or None
        
        :param cut: enforcing the box dimensions are larger than 2x of the cut-off, if not *None*.
        :type cut: float or None
        """
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
                    remained = (now - self.timer_start) / self.reader.nread * (self.reader.frames - self.reader.nread)
                    msg = " -> Processed %d of %d frames. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (self.reader.nread, self.reader.frames, elapsed, remained)
                else:
                    msg = " -> Processed %d frames. Elapsed: %0.0f secs ..." % (self.reader.nread, elapsed)

                screen.info('\r%s' % (msg), end="")

            return self.reader
        else:
            raise StopIteration
        
            
        
        
        
        
        
        
        
        
        
    
