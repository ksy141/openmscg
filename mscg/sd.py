from mscg import *

class SD:
    
    def __init__(self, stepsize = 1.0, maxstep = 1000, ftol = 1e-4):
        self.stepsize, self.maxstep, self.ftol = stepsize, maxstep, ftol
        self.force, self.plist, self.blist = None, None, None
        
    def run(self, X, box):
        F, E = np.zeros(X.shape), 0.0
        prevE = None
        
        for i in range(self.maxstep):
            
            if self.plist is not None:
                self.plist.build(box, X)

            if self.blist is not None:
                self.blist.build(box, X)
        
            F = self.force.compute()
            E = self.force.e
            
            SumF = np.sqrt(np.square(F).sum())
            #print("Step: %d E: %e SumF: %f Step: %f" % (i, E, SumF, self.stepsize))
            
            if prevE is not None and E>prevE:
                self.stepsize *= 0.5
            
            if SumF < self.ftol:
                break
            
            X += F * self.stepsize
            X = Trajectory.pbc(box, X)
            prevE = E
        
        screen.info("SD Minimization Result: step=%d E=%e Sum(F)=%e" % (i, E, SumF))
        return X