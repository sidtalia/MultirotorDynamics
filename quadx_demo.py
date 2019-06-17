'''

  Simple demo of multirotor dynamics Quad-X layout with ArduPilot
  motor-numbering convention:

    3cw   1ccw
       \ /
        X
       / \
    2ccw  4cw

  Copyright (C) 2019 Simon D. Levy

  MIT License

'''

from multirotor_dynamics import Dynamics

class QuadXAP(Dynamics):

    def __init__(self, b, d, m, l, Ix, Iy, Iz, Jr, maxrpm):

        Dynamics.__init__(self, 4, b, d, m, l, Ix, Iy, Iz, Jr, maxrpm)

if __name__ == '__main__':


    # Amir's calculations
    b  = 5.30216718361085E-05,   
    d  = 2.23656692806239E-06,   
    m  = 16.47,                  
    l  = 0.6,                    
    Ix = 2,                      
    Iy = 2,                      
    Iz = 3,                      
    Jr = 3.08013E-04,            

    # estimated
    maxrpm = 15000

    quad = QuadXAP(b, d, m, l, Ix, Iy, Iz, Jr, maxrpm)
