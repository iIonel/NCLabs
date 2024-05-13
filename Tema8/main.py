import random

epsilon = 10 ** -16
niu = 10 ** -5
'''def set_precision_masina():
    t = int(input("Introduceți pentru care 10^(-t) este precizia mașinii: "))
    return t'''

def verificare_epsilon(numar):
    if abs(numar) <= epsilon:
        print(f"Nu se poate împărți, numărul = {numar} este mai mic decât epsilon = {epsilon}!")
        return False
    return True


def prima_derivata(f, x, y, h=1e-6):
    G1 = (3*f(x + h, y) - 4*f(x - h, y)+f(x-2*h,y)) / (2 * h)
    G2 = (3*f(x, y + h) - 4*f(x, y - h)+f(x,y-2*h)) / (2 * h)
    return G1, G2

def a_doua_derivata(f, x, y, h=1e-6):
    d2f_dx2 = (f(x + h, y) - 2 * f(x, y) + f(x - h, y)) / h ** 2
    d2f_dy2 = (f(x, y + h) - 2 * f(x, y) + f(x, y - h)) / h ** 2
    return d2f_dx2, d2f_dy2

def modul_derivata(F, x, y):
    return (prima_derivata(F, x, y)[0] ** 2 + prima_derivata(F, x, y)[1] ** 2)
def calcul_niu(F,x,y):


    niu = 1
    p = 1
    beta = 0.8
    x_derivata = prima_derivata(F, x, y)[0]
    y_derivata = prima_derivata(F, x, y)[1]
    while (F(x-x_derivata,y-y_derivata) > F(x,y) - niu/2 * modul_derivata(F,x,y)) and p < 8:
        niu = beta * niu
        p += 1

    return niu


def initializare_date_random():
    F = lambda x,y: x**2 + y**2 -2*x -4*y -1
    #F = lambda x, y: x**2*y - 2*x*y**2 + 3*x*y + 4
    #F = lambda x, y: x**2 - 4*x*y + 5**y - 4*y + 3
    x0 = random.randint(1, 100)
    y0 = random.randint(1, 100)
    return F, x0, y0

def main():
    k = 0
    F, x0, y0 = initializare_date_random()
    #niu = calcul_niu(F, x0, y0)
    niu = 10
    x = x0 - niu * prima_derivata(F, x0, y0)[0]
    y = y0 - niu * prima_derivata(F, x0, y0)[1]
    k = k + 1
    print(F)
    while verificare_epsilon(niu * modul_derivata(F, x, y)) and (k < 30000) and (niu * modul_derivata(F, x, y) < 10**10):
        x = x - niu * prima_derivata(F, x, y)[0]
        y = y - niu * prima_derivata(F, x, y)[1]
        #niu = calcul_niu(F, x, y)
        niu = 10
        k = k + 1
        print(k)
    if (niu * modul_derivata(F, x, y) > epsilon):
        #print(f"Nu se poate împărți, numărul = {niu * modul_derivata(F, x, y)} este mai mare decât epsilon = {epsilon}!")
        return
    #else:

    print("Solutia: x este", x, "y este", y)

main()