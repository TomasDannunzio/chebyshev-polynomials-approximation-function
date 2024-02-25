import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import scipy as scp

#NOTA: Para que el codigo funcione se debe guardar el archivo dado por la catedra en esta direccion o cambiar la misma.
f = open("C:\Users\tomid\PycharmProjects\TP3Chebyshev", "r")

texto = f.read()

lista = texto.split("\n")

arregloMedidos = []

for i in range(1,len(lista),1):
    x = lista[i].partition(";")[2]
    arregloMedidos.append(float(x))

arregloMedidos.pop(0)

tiempo = np.arange(0, 6, 0.03333)

tiempoLista = tiempo.tolist()

tiempoLista.pop()

print(len(tiempoLista))
print(len(arregloMedidos))

#plt.plot(tiempoLista, arregloMedidos)
plt.xlabel("Tiempo - microsegundos")
plt.ylabel("Potencial Medido - mV")
#plt.show()

z = np.arange(-1, 1, 1/90)

x = sp.symbols('x')

j = []

for i in range(len(z)):
    j.append((1/2)*(180*0.03333)*z[i]+(1/2)*6)

#Funcion que genera todos los polinomios de Chebyshev
def calcularT(orden):
    t0 = sp.sympify(1.)
    t1 = sp.sympify(x)
    funciones = [t0, t1]
    for i in range(orden):
        funciones.append(sp.simplify(sp.sympify(2*x*funciones[1+i]-funciones[i])))
    return funciones

#Este parametro determina el numero de polinomios con los que se aproximaran los datos
n=25

vectorSolucion = []

funciones = calcularT(n)

for i in range(n):
    funcion = []
    for j in range(len(arregloMedidos)):
        funcion.append(funciones[i].subs(x, z[j]))
    vectorSolucion.append(np.sum(arregloMedidos*np.array(funcion)))

print(vectorSolucion)

matrizPrincipal = []

for i in range(n):
    fila = []
    funcioni = []
    for j in range(len(arregloMedidos)):
        funcioni.append(funciones[i].subs(x, z[j]))
    for h in range(n):
        funcionj = []
        for l in range(len(arregloMedidos)):
            funcionj.append(funciones[h].subs(x, z[l]))
        fila.append(np.sum(funcioni*np.array(funcionj)))
    matrizPrincipal.append(fila)

MP = sp.Matrix(matrizPrincipal)
VS = sp.Matrix(vectorSolucion)

solucion = sp.linsolve((MP, VS))

fa = []
for i in range(len(arregloMedidos)):
    fa.append(0)

funcion_aproximante = sp.Matrix(fa)

for i in range(n):
    f = []
    for j in range(len(tiempoLista)):
        f.append(funciones[i].subs(x, (2*tiempoLista[j]-6)/6))
    funcion = sp.Matrix(f)
    funcion_aproximante = funcion_aproximante + solucion.args[0][i] * funcion

print(scp.integrate.simpson(funcion_aproximante[0:68:1], tiempoLista[0:68:1]))

plt.plot(tiempoLista, arregloMedidos)
plt.plot(tiempoLista, funcion_aproximante)
plt.show()