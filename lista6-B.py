from math import *

def f1(V,i):
    x,y,z = V[0][0], V[1][0], V[2][0]
    if (i==1):
        resultado = 16.0*(x**4) + 16.0*(y**4) + 1.0*(z**4) - 16.0
        return resultado
    if (i==2):
        resultado = 1.0*(x**2) + 1.0*(y**2) + 1.0*(z**2) - 3.0
        return resultado
    if (i==3):
        resultado = 1.0*(x**3) - (1.0*y) + 1.0*z - 1.0
        return resultado

def df1(V,i):
    x,y,z = V[0][0], V[1][0], V[2][0]
    if (i==1):
        res1 = 64.0*(x**3)
        res2 = 64.0*(y**3)
        res3 = 4.0*(z**3)
        return [res1, res2, res3]
    if (i==2):
        res1 = 2.0*(x)
        res2 = 2.0*(y)
        res3 = 2.0*z
        return [res1, res2, res3]
    if (i==3):
        res1 = 3.0*(x**2)
        res2 = -1.0
        res3 = 1.0
        return [res1, res2, res3]

def f2(V,i,teta1=0,teta2=11.667):
    c2, c3, c4 = V[0][0], V[1][0], V[2][0]
    if (i==1):
	    ret = 2*(c3**2) + c2**2 + 6*(c4**2) - 1
	    return ret
    if (i==2):
            ret = 8*(c3**3) + 6*c3*(c2**2) + 36*c3*c2*c4 + 108*c3*(c4**2) - teta1
            return ret
    if (i==3):
            ret = 60*(c4**4) + 60*(c3**2)*(c2**2) + 576*(c3**2)*c2*c4 + 2232*(c3**2)*(c4**2) + 252*(c4**2)*(c2**2) + 1296*(c4**3)*c2 + 3348*(c4**4) + 24*(c2**3)*c4 + 3*c2 - teta2
            return ret

def df2(V,i,teta1=0,teta2=11.667):
    c2, c3, c4 = V[0][0], V[1][0], V[2][0]
    #derivadas em relacao a c2
    if (i==1):
        res1 = 2*c2
        res2 = 4*c3
        res3 = 12*c4
        return [res1, res2, res3]
    #derivadas em relacao a c3
    if (i==2):
        res1 = 12*c3*c2 + 36*c3*c4
        res2 = 24*(c3**2) + 6*(c2**2) + 108*(c4**2)
        res3 = 36*c3*c2 + 216*c3*c4
        return [res1, res2, res3]
    #derivadas em relacao a c4
    if (i==3):
        res1 = 20*(c3**2)*c2 + 576*(c3**2)*c4 + 504*(c4**2)*c2 + 1296*(c4**3)
        + 72*(c2**2)*c4 + 3
        res2 = 240*(c3**3) + 120*c3*(c2**2) + 1152*c3*c2*c4 + 4464*c3*(c4**2)
        res3 = 576*(c3**2)*c2 + 4464*(c3**2)*c4 + 504*c4*(c2**2) + 3888*(c4**2)*c2 + 13392*(c4**3) + 24*(c2**3)
        return [res1, res2, res3]

def f3(V,i):
    b0,b1,b2 = V[0][0], V[1][0], V[2][0]
    x = [1.0,2.0,3.0]
    y = [1.0,2.0,9.0]
    ret = b0 + b1*(x[i-1]**b2) - y[i-1]
    return ret

def df3(V,i):
    b0,b1,b2 = V[0][0], V[1][0], V[2][0]
    x = [1.0,2.0,3.0]
    ret1 = 1.0
    ret2 = x[i-1]**b2
    ret3 = (b1*(x[i-1]**b2)*log(x[i-1]))
    return [ret1, ret2, ret3]

def funcao(nfuncao, V, i):
    if nfuncao<1 or nfuncao>3:
        return "Numero da funcao invalido"
    if (nfuncao==1):
        return f1(V, i)
    if (nfuncao==2):
        return f2(V, i)
    if (nfuncao==3):
        return f3(V, i)

def dfuncao(nfuncao, V, i):
    if nfuncao<1 or nfuncao>3:
        return "Numero da funcao invalido"
    if (nfuncao==1):
        return df1(V, i)
    if (nfuncao==2):
        return df2(V, i)
    if (nfuncao==3):
        return df3(V, i)

def euclidean(X):
    total = 0;
    for x in X:
        total += 1.0*x[0]**2
    return sqrt(total)
                  
def transpose(M):
    return [list(i) for i in zip(*M)] #unzipping para uma lista

def getMatrixMinor(M,i,j):
    ret = [row[:j] + row[j+1:] for row in (M[:i]+M[i+1:])]
    return ret

def getDeterminant(M):
    if len(M)==1:
        return M[0][0]
    #caso especial para matrizes 2x2
    if len(M) == 2:
        return M[0][0]*M[1][1]-M[0][1]*M[1][0]
    #
    determinant = 0
    for i in range(len(M)):
        determinant += ((-1)**i)*M[0][i]*getDeterminant(getMatrixMinor(M,0,i))
    return determinant

def getInverse(M):
    determinant = getDeterminant(M)
    if determinant==0:
        print "Matrix cannot be inversed. Its determinant is 0."
        return "Matrix cannot be inversed. Its determinant is 0."
    if len(M)==1:
        return [[1/M[0][0]]]
    #caso especial para matrizes 2x2
    if len(M) == 2:
        return [[M[1][1]/determinant, -1*M[0][1]/determinant],
                [-1*M[1][0]/determinant, M[0][0]/determinant]]
    #encontrar matriz de cofatores
    cofactors = []
    for i in range(len(M)):
        cofactorRow = []
        for j in range(len(M)):
            minor = getMatrixMinor(M,i,j)
            cofactorRow+=[((-1)**(i+j)) * getDeterminant(minor)]
        cofactors+=[cofactorRow]
    cofactors = transpose(cofactors)
    for i in range(len(cofactors)):
        for j in range(len(cofactors)):
            cofactors[i][j] = cofactors[i][j]/determinant
    return cofactors

def mult(A,B):
    #height: qtd de linhas
    #width: qtd de colunas
    heightA, widthA = len(A), len(A[0])
    heightB, widthB = len(B), len(B[0])
    #verificar dimensoes das matrizes
    if (widthA != heightB):
        return "Dimensoes das matrizes sao incompativeis"
    #inicializar matriz Resultado
    R = [[0.0 for x in range(widthB)] for y in range(heightA)] 
    #
    for i in range(heightA):
        for j in range(widthB):
            for k in range(widthA):
                R[i][j] += A[i][k] * B[k][j]
    return R

def add(A,B):
    R = []
    for i in range(len(A)):
        row = []
        for j in range(len(A[i])):
            row+=[A[i][j] + B[i][j]]
        R+=[row]
    return R

def scalar(a,A):
    ret = [[(a)* x for x in A[i]] for i in range(len(A))]
    return ret

def NewtonMethod(nfuncao, X0=[[1.0],[1.0],[1.0]], n=3, niter=1000, tol=10.0**(-4)): #X0=[[-1.0],[0.0],[1.0]]
    tolk = 0.0
    for k in range (niter):
        J, F, deltaX = [], [], []
        for i in range(n):
            F += [[funcao(nfuncao,X0,i+1)]]
            J += [dfuncao(nfuncao,X0,i+1)]
        Ji = getInverse(J)
        deltaX = scalar(-1,mult(Ji,F))
        X0 = add(X0,deltaX)
        tolk = 1.0*euclidean(deltaX)/euclidean(X0)
        if (tolk < tol):
            return X0
    return "Convergencia NAO ALCANCADA"

def BroydenMethod(nfuncao, X0=[[1.0],[1.0],[1.0]], B=[[0,-3,-1],[3,0,5],[1,5,0]], niter=1000, tol=10.0**(-4)):
    tolk = 0.0
    n=len(B)
    for k in range (niter):
        Y, F, deltaX = [], [], []
        J = B
        for i in range(n):
            F += [[funcao(nfuncao,X0,i+1)]]
        #inverter a matriz para resolver o sistema
        Ji = getInverse(J)
        #multiplicar matriz inversa Ji pelo vetor V0
        deltaX = scalar(-1,mult(Ji,F))
        #salvar valor antigo de F
        F0 = F
        #atualizar X0
        X0 = add(X0,deltaX)
        #calcular F com X0 atualizado
        F = []
        for i in range(n):
            F += [[funcao(nfuncao,X0,i+1)]]
        #calcular Y
        Y = add(F,scalar(-1,F0))
        #calcular o erro e comparar com tolerancia
        tolk = 1.0*euclidean(deltaX)/euclidean(X0)
        if (tolk < tol):
            return X0
        else:
            aux = mult(B,deltaX)
            aux = add(Y,scalar(-1.0,aux))
            aux = mult(aux,transpose(deltaX))
            aux2 = mult(transpose(deltaX), deltaX)
            B = add(B,scalar(getInverse(aux2)[0][0], aux))
    return "Convergencia NAO ALCANCADA"

def NL_MinimumSquares(nfuncao, B0=[[3.0],[2.0],[1.0]], n=3, niter=1000, tol=10.0**(-4)):
    tolk = 0.0
    #B0 = [[1.0],[2.0],[3.0]]
    for k in range (niter):
        J, F, deltaB = [], [], []
        for i in range(n):
            F += [[funcao(nfuncao,B0,i+1)]]
            J += [dfuncao(nfuncao,B0,i+1)]
        #calcular deltaB: -inv(J'*J)*J'*F
        a = mult(transpose(J),J)
        a = getInverse(a)
        b = mult(transpose(J),F)
        deltaB = scalar(-1,mult(a,b))
        #Calcular B0 como no metodo de Newton
        B0 = add(B0,deltaB)
        tolk = 1.0*euclidean(deltaB)/euclidean(B0)
        if (tolk < tol):
            return B0
    return "Convergencia NAO ALCANCADA"
