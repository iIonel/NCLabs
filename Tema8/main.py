import random

epsilon = 10 ** -16
niu = 10 ** -5

def verificare_epsilon(numar):
    if abs(numar) <= epsilon:
        print(f"Nu se poate împărți, numărul = {numar} este mai mic decât epsilon = {epsilon}!")
        return False
    return True

def prima_derivata(f, x, y, h=1e-6):
    derivata_x = (3*f(x + h, y) - 4*f(x - h, y) + f(x - 2*h, y)) / (2 * h)
    derivata_y = (3*f(x, y + h) - 4*f(x, y - h) + f(x, y - 2*h)) / (2 * h)
    return derivata_x, derivata_y

def a_doua_derivata(f, x, y, h=1e-6):
    derivata2_x = (f(x + h, y) - 2 * f(x, y) + f(x - h, y)) / h ** 2
    derivata2_y = (f(x, y + h) - 2 * f(x, y) + f(x, y - h)) / h ** 2
    return derivata2_x, derivata2_y

def modul_derivatei(F, x, y):
    derivata_x, derivata_y = prima_derivata(F, x, y)
    return derivata_x ** 2 + derivata_y ** 2

def calcul_niu(F, x, y):
    niu = 1
    p = 1
    beta = 0.8
    derivata_x, derivata_y = prima_derivata(F, x, y)
    while (F(x - derivata_x, y - derivata_y) > F(x, y) - niu / 2 * modul_derivatei(F, x, y)) and p < 8:
        niu = beta * niu
        p += 1
    return niu

def initializare_date_random():
    F = lambda x, y: x**2 + y**2 - 2*x - 4*y - 1
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
    derivata_x, derivata_y = prima_derivata(F, x0, y0)
    x = x0 - niu * derivata_x
    y = y0 - niu * derivata_y
    k = k + 1
    print(F)
    while verificare_epsilon(niu * modul_derivatei(F, x, y)) and (k < 30000) and (niu * modul_derivatei(F, x, y) < 10**10):
        derivata_x, derivata_y = prima_derivata(F, x, y)
        x = x - niu * derivata_x
        y = y - niu * derivata_y
        #niu = calcul_niu(F, x, y)
        niu = 10
        k = k + 1
        print(k)
    if niu * modul_derivatei(F, x, y) > epsilon:
        return
    print("Soluția: x este", x, "y este", y)

main()