from math import *

def bissecao(nFuncao, a, b, tol=0.0001):
    """Lista 4: (1,250,300) (1,-300,-250) (2,-10,-5) (2,0,2)"""
    if a==b:
        if funcao(nFuncao, a-tol)<0<funcao(nFuncao, a+tol):
            return a
        else:
            return "Parametro nao e raiz. Insira dois parametros diferentes"
    if funcao(nFuncao, a)>funcao(nFuncao, b):
        a,b=b,a
    while abs(b-a)>=tol:
        x=(a+b)/2.0
        f=funcao(nFuncao, x)
        if f>0:
            b=x
        else:
            a=x
    return x

def newtonRaiz(nFuncao, x0, secante=False, ITERMAX=1000000, tol=0.0001):
    if not secante:
        for i in range(ITERMAX):
            x=x0-funcao(nFuncao, x0)/dfuncao(nFuncao, x0)
            if abs(x-x0)<tol:
                return x
            x0=x
    else:
        delta=0.001
        x=x0+delta
        fa=funcao(nFuncao, x0)
        for i in range(ITERMAX):
            fi=funcao(nFuncao, x)
            xprox=x-fi*(x-x0)/(fi-fa)
            if abs(xprox-x)<tol:
                return xprox
            fa=fi
            x0=x
            x=xprox
    return "Convergencia nao foi alcancada"

def interpolacaoInversa(nFuncao, X, ITERMAX=1000000, tol=0.0001):
    if len(X)<3:
        return "Insert more points"
    if len(set(X))<len(X):
        return "All points must be distinct"
    X.sort()
    x0=None
    for i in range(ITERMAX):
        Y=[]
        for i in X:
            Y+=[funcao(nFuncao, i)]
        x=0
        for i in range(len(X)):
            numerador=1
            denominador=1
            for k in range(len(Y)):
                if k!=i:
                    numerador*=Y[k]
                    denominador*=Y[i]-Y[k]
            x+=X[i]*numerador/denominador
        if x0!=None and abs(x-x0)<tol:
            return x
        x0=x
        if abs(max(Y))>abs(min(Y)):
            i=Y.index(max(Y))
        else:
            i=Y.index(min(Y))
        X[i]=x
        X.sort()
    return "Convergencia nao foi alcancada"

def funcao(n, x):
    if n<1 or n>2:
        return "Funcao desconhecida"
    if n==1:
        g=9.806
        k=0.00341
        return log(cosh(x*sqrt(g*k)))-50
    if n==2:
        return 4*cos(x)-exp(2*x)

def dfuncao(n, x):
    if n<1 or n>2:
        return "Funcao desconhecida"
    if n==1:
        g=9.806
        k=0.00341
        return sqrt(g*k)*tanh(x*sqrt(g*k))
    if n==2:
        return -2*(exp(2*x)+2*sin(x))

def funcao1(x, g=9.806, k=0.00341):
    return log(cosh(x*sqrt(g*k)))-50

def funcao2(x):
    return 4*cos(x)-exp(2*x)

def dfuncao1(x, g=9.806, k=0.00341):
    return sqrt(g*k)*tanh(x*sqrt(g*k))

def dfuncao2(x):
    return -2*(exp(2*x)+2*sin(x))
