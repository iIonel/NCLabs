import math
import random


class Muller:
    def __init__(self, coeficienti_polinom):
        self.coeficienti_polinom = coeficienti_polinom

    def horner(self, coeficienti, x):
        rezultat = coeficienti[0]
        for coef in coeficienti[1:]:
            rezultat = rezultat * x + coef
        return rezultat

    def calculeaza_abc(self, coeficienti, x):
        x_0, x_1, x_2 = x
        h_0 = x_1 - x_0
        h_1 = x_2 - x_1
        delta_0 = self.calculeaza_delta(coeficienti, h_0, x_1, x_0)
        delta_1 = self.calculeaza_delta(coeficienti, h_1, x_2, x_1)
        a = self.calculeaza_div((delta_1 - delta_0), (h_1 + h_0))
        b = a * h_1 + delta_1
        c = self.horner(coeficienti, x_2)
        return {"x_0": x_0, "x_1": x_1, "x_2": x_2, "h_0": h_0, "h_1": h_1, "delta_0": delta_0, "delta_1": delta_1,
                "a": a, "b": b, "c": c}

    def calculeaza_polinom(self, coeficienti, x):
        rezultat = coeficienti[0]
        for coef in coeficienti[1:]:
            rezultat = rezultat * x + coef
        return rezultat

    def calculeaza_delta(self, coeficienti, h, x, y):
        return self.calculeaza_div((self.calculeaza_polinom(coeficienti, x) - self.calculeaza_polinom(coeficienti, y)),
                                   h)

    def calculeaza_div(self, a, b):
        if b == 0:
            return 0
        return a / b

    def semn(self, x):
        if x >= 0:
            return 1
        else:
            return -1

    def metoda_muller(self, f, x0, x1, x2, tol=1e-6, max_iter=100):
        numar_iteratii = 0
        while numar_iteratii < max_iter:
            numar_iteratii += 1
            rezultat = self.calculeaza_abc(self.coeficienti_polinom, [x0, x1, x2])
            delta = rezultat['b'] ** 2 - 4 * rezultat['a'] * rezultat['c']
            if delta < 0:
                return None
            delta_x = self.calculeaza_div(2 * rezultat['c'],
                                          rezultat['b'] + self.semn(rezultat['b']) * math.sqrt(delta))
            x_3 = x2 - delta_x
            if abs(delta_x) < tol:
                return x_3
            x0, x1, x2 = x1, x2, x_3
        return None

    def gaseste_radacini(self, puncte_start, tol=1e-6, max_iter=100):
        radacini = set()
        for i in range(len(puncte_start)):
            x0, x1, x2 = puncte_start[i], puncte_start[(i + 1) % len(puncte_start)], puncte_start[
                (i + 2) % len(puncte_start)]
            radacina = self.metoda_muller(lambda x: self.horner(self.coeficienti_polinom, x), x0, x1, x2, tol, max_iter)
            if radacina is not None:
                radacini.add(round(radacina, 6))
        return radacini

    def calculeaza_interval_radacini(self):
        A = max(abs(coef) for coef in self.coeficienti_polinom[1:])
        R = (abs(self.coeficienti_polinom[0]) + A) / abs(self.coeficienti_polinom[0])
        return -R, R


if __name__ == "__main__":
    coeficienti_polinom = [1.0, -6.0, 11.0, -6.0]
    muller = Muller(coeficienti_polinom)

    interval_radacini = muller.calculeaza_interval_radacini()
    print("Interval radacini: ", interval_radacini)

    radacini = set()
    for i in range(100):
        puncte_start = [random.uniform(interval_radacini[0], interval_radacini[1]) for _ in range(3)]
        radacini.update(muller.gaseste_radacini(puncte_start))

    with open("radacini.txt", "w") as fisier:
        for radacina in radacini:
            fisier.write(str(radacina) + "\n")