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

    # roll right
    def u2(o):
    
        return (o[1] + o[2]) - (o[0] + o[3])
    

    # pitch forward
    def u3(o):
    
        return (o[1] + o[3]) - (o[0] + o[2])
    

    # yaw cw
    def u4(o):
    
        return (o[0] + o[1]) - (o[2] + o[3])
