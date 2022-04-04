def response(t, f):

    #Initial Conditions:
    x0=0 
    v0=0

    nPoints = len(t)  
    n = 0
    pos = [] 
    vel = []

    for timestep in t:
        
        if (n<len(t)-1):
            Fn = f[n] 
            Fnp1 = f[n+1]
        else:
            Fn = f[n]
            Fnp1 = 0

        if (n<len(t)-1):
            curPos = (A*x0) + (B*v0) + (C*Fn) + (D*Fnp1) 
            curVel = (A1*x0) + (B1*v0)+ (C1*Fn) + (D1*Fnp1)
            
            pos.append(curPos) 
            vel.append(curVel) 

            #Update Initial Conditions:
            x0 = curPos;
            v0 = curVel;

            n = n + 1

        else:
            curPos = (A*x0) + (B*v0) + (C*Fn) + (D*0) 
            curVel = (A1*x0) + (B1*v0)+ (C1*Fn) + (D1*0)
            
            pos.append(curPos) 
            vel.append(curVel) 

            #Update Initial Conditions:
            x0 = curPos;
            v0 = curVel;

            n = n + 1

    return pos