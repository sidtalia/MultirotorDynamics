'''
  Dynamics for a big (16kg) quadcopter

  Copyright (C) 2019 Simon D. Levy

  MIT License

'''

from quadxap import QuadXAP

class BigQuad(QuadXAP):

    # Amir's calculations
    b  = 5.30216718361085E-05   
    d  = 2.23656692806239E-06   
    m  = 16.47                  
    l  = 0.6                    
    Ix = 2                      
    Iy = 2                      
    Iz = 3                      
    Jr = 3.08013E-04            

    # estimated
    maxrpm = 15000

    def __init__(self):

        QuadXAP.__init__(self, self.b, self.d, self.m, self.l, self.Ix, self.Iy, self.Iz, self.Jr, self.maxrpm)
