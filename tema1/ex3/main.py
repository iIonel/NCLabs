import random
import math

t_calcul = [
    lambda a: (105 * a - 10 * (a ** 3)) / (105 - 45 * (a ** 2) + a ** 4), #t4
    lambda a: (945 * a - 105 * (a ** 3) + a ** 5) / (945 - 420 * (a ** 2) + 15 * (a ** 4)), #t5
    lambda a: (10395 - 1260 * (a ** 3) + 21 * (a ** 5)) / (10395 - 4725 * (a ** 2) + 210 * (a ** 4) - (a ** 6)), #t6
    lambda a: (135135 * a - 17325 * (a ** 3) + 378 * (a ** 5) - (a ** 7)) / (135135 - 62370 * (a ** 2) + 3150 * (a ** 4) - 28 * (a ** 6)), #t7
    lambda a: (2027025 * a - 270270 * (a ** 3) + 6930 * (a ** 5) - 36 * (a ** 7)) / (2027025 - 945945 * (a ** 2) + 51975 * (a ** 4) - 630 * (a ** 6) + (a ** 8)), #t8
    lambda a: (34459425 * a - 4729725 * (a ** 3) + 135135 * (a ** 5) - 990 * (a ** 7) + (a ** 9)) / (34459425 - 16216200 * (a ** 2) + 945945 * (a ** 4) - 13860 * (a ** 6) + 45 * (a ** 8))] #t9

def gen_numbers():
    numbers = []
    for i in range(1, 10000 + 1):
        x = random.uniform(-math.pi / 2, math.pi / 2)
        numbers.append(x)
    return numbers

def aprox():
    with open("results.txt", "w") as file:
        counting = {4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}

        test_numbers = gen_numbers()
        file.write(f"Random numbers [-pi/2, pi/2]: {test_numbers}")
        file.write("\n")

        for number in test_numbers:
            t_aprox = [(i + 4,float(abs(float(t_calcul[i](number)) - float(math.tan(number))))) for i in range(0,6)]
            top = sorted(t_aprox,key=lambda x:x[1])
            position = 0
            for t in top:
                counting[t[0]] += position
                position = position + 1
            file.write(f"Top3: {top[:3]}")
            file.write("\n")

        final = [(k,v) for k,v in counting.items()]
        file.write("\n")
        file.write(f"{final}")
        final.sort(key = lambda x: x[1])
        final = [k for (k,v) in final]
        file.write("\n")
        file.write(f"Ranking: {final}")

        file.write("\n")
        file.write("\n")

        for number in test_numbers:
            file.write(f"Number {number}: T6 = {t_calcul[2](number)}, S6 = {t_calcul[2](number)/math.sqrt(1 + (t_calcul[2](number)) * (t_calcul[2](number)))}, C6 = {1/math.sqrt(1 + (t_calcul[2](number)) * (t_calcul[2](number)))}, T7 = {t_calcul[3](number)}, S7 = {t_calcul[3](number)/math.sqrt(1 + (t_calcul[3](number)) * (t_calcul[3](number)))}, C7 = {1/math.sqrt(1 + (t_calcul[3](number)) * (t_calcul[3](number)))} ")
            file.write("\n")

if __name__ == '__main__':
    aprox()
