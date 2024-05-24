m = 1
u = pow(10, -m)

def exercice1():
    global m,u

    while 1 + u != 1:
        m += 1
        u = pow(10, -m)
    m = m - 1

    print(m)
    print (u)

    return u
    
exercice1()
