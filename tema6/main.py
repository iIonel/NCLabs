import numpy as np
import matplotlib.pyplot as plt

def interpolare_aitken(valori_x, valori_y, x_tinta):
    numar_puncte = len(valori_x)
    matrice_F = np.zeros((numar_puncte, numar_puncte))
    for i in range(numar_puncte):
        matrice_F[i, 0] = valori_y[i]
    for j in range(1, numar_puncte):
        for i in range(numar_puncte - j):
            matrice_F[i, j] = matrice_F[i + 1, j - 1] + (matrice_F[i + 1, j - 1] - matrice_F[i, j - 1]) / (valori_x[i + j] - valori_x[i])
    for k in range(numar_puncte):
        if valori_x[k] == x_tinta:
            return matrice_F[0, k]
    return None

def functie_tinta(x):
    return x ** 4 - 12 * x ** 3 + 30 * x ** 2 + 12

def aproximare_patratica(valori_x, valori_y, x_tinta):
    numar_puncte = len(valori_x)
    for grad in range(1, min(6, numar_puncte)):
        matrice_vandermonde = np.vander(valori_x, grad + 1, increasing=True)
        coeficienti, _, _, _ = np.linalg.lstsq(matrice_vandermonde, valori_y, rcond=None)
        valoare_polinom = np.polyval(coeficienti, x_tinta)
        return valoare_polinom

if __name__ == '__main__':
    numar_puncte = int(input("Introdu numărul de puncte: "))
    x_initial = float(input("x0: "))
    x_final = float(input("xn: "))
    pas = (x_final - x_initial) / numar_puncte

    valori_x = [x_initial + i * pas for i in range(numar_puncte + 1)]
    valori_y = [functie_tinta(x) for x in valori_x]
    print("Valori x: ", valori_x)
    print("Valori y: ", valori_y)

    x_nou = float(input("Introdu valoarea nouă a lui x: "))

    rezultat_aitken = interpolare_aitken(valori_x, valori_y, x_nou)
    rezultat_patratica = aproximare_patratica(valori_x, valori_y, x_nou)

    print("Rezultat Interpolare Aitken: ", rezultat_aitken)
    print("Rezultat Aproximare Pătratică: ", rezultat_patratica)

    if rezultat_aitken is not None:
        print("|Ln(x) - f(x)|: ", abs(rezultat_aitken - functie_tinta(x_nou)))

    plt.plot(valori_x, valori_y, label='f(x)', color='black', linestyle='dashed')
    if rezultat_aitken is not None:
        plt.plot(x_nou, rezultat_aitken, 'bo', label='Ln(x)')
    plt.plot(x_nou, rezultat_patratica, 'ro', label='Pm(x)')

    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()
