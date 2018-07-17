#!/usr/bin/python3
# -*- coding: UTF-8 -*-

def blue(J,D,tRL,tLR,UR,UL,e1,e2):
   """
   order: {u,u},{d,d},{u,d},{d,u},{ud,},{,ud}
   """
   M = np.matrix([[-J/4,  0 ,  0 ,  0 ,   0  ,   0  ],
                  [  0 ,-J/4,  0 ,  0 ,   0  ,   0  ],
                  [  0 ,  0 , J/4,-J/2,  tRL ,  tLR ],
                  [  0 ,  0 ,-J/2, J/4, -tRL , -tLR ],
                  [  0 ,  0 , tRL,-tRL,UR-J/4,   D  ],
                  [  0 ,  0 , tLR,-tLR,  D   ,UL-J/4]])

   K = np.matrix([[e1+e2,  0  ,  0  ,  0  ,  0  ,  0  ],
                  [  0  ,e1+e2,  0  ,  0  ,  0  ,  0  ],
                  [  0  ,  0  ,e1+e2,  0  ,  0  ,  0  ],
                  [  0  ,  0  ,  0  ,e1+e2,  0  ,  0  ],
                  [  0  ,  0  ,  0  ,  0  ,e1+e1,  0  ],
                  [  0  ,  0  ,  0  ,  0  ,  0  ,e2+e2]])
   return M+K
