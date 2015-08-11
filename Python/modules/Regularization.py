from scipy import optimize
import numpy as np
import math
from pyglow import pyglow

def get_regularized_ed(lamda,Rayl2=0.,NE=0.,S=0.,rbot=0.):
    
    
    # Construct the Diagonal Identity matrix containing the regularization parameter 
    x = np.zeros((np.size(S,0),np.size(S,0)))
    x[:] = lamda
    I = np.diag(np.diag(x))

    # Constuct matrix C= [S;lambda*I]
    C = np.concatenate((S, I))

    # Constuct vector d = [b;0]
    b = np.zeros(np.size(Rayl2))
    d = np.concatenate((Rayl2, b))

    sol = optimize.nnls(C,d)

    res = sol[0]
    a1356 = 7.3e-13
    counts = np.sqrt(res/a1356)

    # First Order
    L1 = np.zeros((np.size(S,0),np.size(S,0)))
    for i in range(0,(np.size(S,0))):
        L1[i][i] = -1
        if i<(np.size(S,0)-1):
            L1[i][i+1] = 1

    # Constuct matrix C= [S;lambda*I]
    C = np.concatenate((S, lamda*L1))

    # Constuct vector d = [b;0]
    b = np.zeros(np.size(Rayl2))
    d = np.concatenate((Rayl2, b))

    sol = optimize.nnls(C,d)

    res = sol[0]
    a1356 = 7.3e-13
    counts_L1 = np.sqrt(res/a1356)

    # Second Order
    L2 = np.zeros((np.size(S,0),np.size(S,0)))
    for i in range(0,(np.size(S,0))):
        L2[i][i] = 1
        if i<(np.size(S,0)-1):
            L2[i][i+1] = -2
        if i<(np.size(S,0)-2):
            L2[i][i+2] = 1

    # Constuct matrix C= [S;lambda*I]
    C = np.concatenate((S, lamda*L2))

    # Constuct vector d = [b;0]
    b = np.zeros(np.size(Rayl2))
    d = np.concatenate((Rayl2, b))

    sol = optimize.nnls(C,d)

    res = sol[0]
    a1356 = 7.3e-13
    counts_L2 = np.sqrt(res/a1356) 
    
    return counts, counts_L1, counts_L2

def guassian_elinimation(S,Rayl):
    Ip = np.zeros(np.size(Rayl))
    Ip = np.linalg.solve(S,Rayl)
    
    # Remove any negative values
    for i in range(0,np.size(Ip)):
        if Ip[i]<0:
            Ip[i] = 0
    return Ip

def calc_electron_density(Ip,S,O,cont=1):
    '''
    cont = 2 => only RR
    cont = 1 => RR+MN
    '''
    if np.size(np.shape(Ip))==2:
        Ip = Ip[:,0]
    elif np.size(np.shape(Ip))>2:
        raise Exception('Invalid data vector')

    b1356 = 0.54    # yield parameter (unitless)
    a1356 = 7.3e-13 # radiative recombination rate (cm^3/s)
    # consider that to be constant through altitude, normally it should change with temperature
    k1 = 1.3e-15    # radiative attachment rate (cm^3/s)
    k2 = 1e-7       # ion-ion neutralization rate (cm^3/s)
    k3 = 1.4e-10    # ion-atom neutralization rate (cm^3/2)
    
    if cont==1:
        a0 = a1356/O
        b0 = a1356*k3/k2+b1356*k1
        c0 = -Ip/O
        d0 = -Ip*k3/k2

        a1 = b0*c0/(6.0*a0**2)-b0**3./a0**3/27.0-d0/a0/2.0
        
        b1 = c0/a0/3.0-b0**2./a0**2/9.0

        c1 = -b0/a0/3.0;
        
        d1 = a1/np.sqrt((-b1)**3)
        
        d1[np.where(d1<-1.)]=-1.
        d1[np.where(d1>1.)]=1.
            
        # used arccos instead of MATLAB acos
        NE_est = c1+2.0*np.sqrt(-b1)*np.cos(np.arccos(d1)/3.0)
    else:    
        NE_est = np.sqrt(Ip/a1356)

    return NE_est

# New Regularization Methods Post ECE558 Code

# Create vector containing regularization parameters
def create_alpha_values(A,npoints = 1000.):
    # SVD Decomposition of matrix A (Distance matrix)
    U, s, V = np.linalg.svd(A, full_matrices=True)
    
    # multiplication ratio
    smin_ratio = 16*np.spacing(1)
    reg_param = np.zeros(npoints)
    reg_param[npoints-1] = max([s[np.size(s,0)-1],s[1]*smin_ratio])
    ratio = (s[0]*100/reg_param[npoints-1])**(1./(npoints-1));
    
    # Put regularization parameters in descending order
    for i in np.arange(npoints-2,-1,-1):
        reg_param[i] = ratio*reg_param[i+1]
        
    return reg_param

# Find optimum reg param using Maximum Curvature of L-Curve
def Maximum_Curvature(residual,x_lamda,reg_param):

    # Maximum Curvature
    #transform rho and eta into log-log space
    x=np.log(residual);
    y=np.log(x_lamda);

    # the series of points used for the triangle/circle
    x1 = x[:-2]
    x2 = x[1:-1]
    x3 = x[2:]
    y1 = y[:-2]
    y2 = y[1:-1]
    y3 = y[:-2]

    # the side lengths for each triangle
    a1 = np.sqrt((x3-x2)**2+(y3-y2)**2);
    b1 = np.sqrt((x1-x3)**2+(y1-y3)**2);
    c1 = np.sqrt((x2-x1)**2+(y2-y1)**2);

    # semi-perimeter
    s1=(a1+b1+c1)/2;

    # the radius of each circle
    R=(a1*b1*c1)/(4*np.sqrt((s1*(s1-a1)*(s1-b1)*(s1-c1))));

    # The curvature for each estimate for each value which is
    # the reciprocal of its circumscribed radius. Since there aren't circles for 
    # the end points they have no curvature
    kappa = np.zeros(np.size(R)+2)
    kappa[0] = 0
    kappa[-1] = 0
    kappa[1:-1] = 1/R

    ireg_corner=np.argmax(abs(kappa[1:-1]));
    reg_corner=reg_param[ireg_corner];

    return reg_corner

# Tikhonov Regularization
def Tikhonov(A,b,deg,reg_param=0.,ireg=0):
    
    # Make a change to include the shape of reg_param so we dont have to calculate it each time
    # Check if reg_param is a vector or not
    if np.size(reg_param) == 1:
        # if its not a vector only if zero it will calculate the vector. Otherwise it will use the single value given. 
        if reg_param == 0:
            reg_param = create_alpha_values(A)
    #if reg_param == 0:
        #reg_param = create_alpha_values(A)

    residual = np.zeros(len(reg_param))
    seminorm = np.zeros(len(reg_param))
    
    L = get_rough_matrix(len(b),deg)
    # This is the big bottleneck on my code
    for i in range(0,len(reg_param)):
        sol = calc_solution(A,b,reg_param[i],L) 
        residual[i] = np.linalg.norm(b - A.dot(sol))
        seminorm[i] = np.linalg.norm(L.dot(sol))
        
        
    reg_corner = Maximum_Curvature(residual,seminorm,reg_param)
    sol = calc_solution(A,b,reg_corner,L) 
    
    print 'Alpha parameter chose:.%4f' %(reg_corner)

    if ireg==0:
        return sol
    else:
        return sol,reg_corner

# Find optimum regularization parameter using GCV
def gcv(A,b,deg,reg_param=0):


    # Make a change to include the shape of reg_param so we dont have to calculate it each time
    # Check if reg_param is a vector or not
    if np.size(reg_param) == 1:
        # if its not a vector only if zero it will calculate the vector. Otherwise it will use the single value given. 
        if reg_param == 0:
            reg_param = create_alpha_values(A)
    #if reg_param == 0:
        #reg_param = create_alpha_values(A)

    m = len(b)
    I = np.eye(np.size(A,0))
    g = np.zeros(len(reg_param))
    
    L = get_rough_matrix(len(b),deg)

    for i in range(0,len(reg_param)):
        A_sharp =  np.linalg.inv(((A.T).dot(A)+reg_param[i]**2*(L.T).dot(L))).dot(A.T)
        '''
        C = np.concatenate((A, reg_param[i]*L))
        
        # Constuct vector d = [b;0]
        d_temp = np.zeros(len(b))
        d = np.concatenate((b, d_temp))
       
        # Might need to add augmentation here
        sol,rnorm = optimize.nnls(C,d)       
        '''
        sol = calc_solution(A,b,reg_param[i],L) 
        residual = np.linalg.norm(b - A.dot(sol))**2 
        #residual[i] = rnorm
        denom = np.trace(I - A.dot(A_sharp))**2

        g[i] = m*residual/denom

    ireg_corner = np.argmin(g)
    reg_corner=reg_param[ireg_corner];

    print reg_corner

    sol = calc_solution(A,b,reg_corner,L) 

    return sol

# Calculate roughening matrices 1st and 2nd degree
def get_rough_matrix(n,deg):
    if deg==0:
        L = np.eye(n)
        return L
    elif deg==1:
            # First Order
        L = np.zeros((n,n))
        for i in range(0,n):
            L[i][i] = -1
            if i<(n-1):
                L[i][i+1] = 1
        return L
    elif deg==2:
        L = np.zeros((n,n))
        for i in range(0,n):
            L[i][i] = 1
            if i<(n-1):
                L[i][i+1] = -2
            if i<(n-2):
                L[i][i+2] = 1
        return L
    else:
        raise Exception('Invalid degree')

# Calculate regularized inverse rolution
def calc_solution(A,b,lamda,L):
    
    C = np.concatenate((A, lamda*L))
    if np.size(np.shape(b))==2:
        b = b[:,0]
    elif np.size(np.shape(b))>2:
        raise Exception('Invalid data vector')
        
    # Constuct vector d = [b;0]
    d_temp = np.zeros(len(b))
    d = np.concatenate((b, d_temp))

    sol,rnorm = optimize.nnls(C,d)
    
    return sol
                


