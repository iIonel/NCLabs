x = 1.0
y = 1e-16
z = 1e-16
def exercice2():
    global x,y,z

    l = (x+y)+z
    r = x+(y+z)

    if l == r:
        print("asociativa")
    else:
        print("neasociativa")

    x = 1.0e1
    l = (x * y) * z
    r = x * (y * z)

    if l == r:
        print("asociativa")
    else:
        print("neasociativa")

exercice2()