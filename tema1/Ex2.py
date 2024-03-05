def Ex2():
    x = 1.0
    y = 1e-16
    z = 1e-16
    Left = (x+y)+z
    Right = x+(y+z)
    if (Left == Right):
        print("Adunarea este asociativa")
    else:
        print("Adunarea nu este asociativa")
    x = 1.0e1
    Left = (x * y) * z
    Right = x * (y * z)
    if (Left == Right):
        print("Inmultirea este asociativa")
    else:
        print("Inmutirea nu este asociativa")

Ex2()
