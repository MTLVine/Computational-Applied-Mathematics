import math
import matplotlib.pyplot as plt

def rCubed(x, y):
    return pow(pow(x, 2) + pow(y, 2), 3 / 2)

def model(e):
    xModel, yModel, h = [], [], math.pi / 500
    for i in range(1000):
        x, y = (math.cos(i * h) + e), math.sqrt(1 - pow(e, 2)) * math.sin(i * h)
        xModel.append(x)
        yModel.append(y)
    return (xModel, yModel)

def FE(d, e, N):
    x, y, u, v = 1 + e, 0, 0, math.sqrt((1 - e) / (1 + e))
    h = 2 * math.pi / d
    xValues, yValues = [x], [y]
    while len(xValues) < d * N:
        xNew, yNew = x + h * u, y + h * v
        uNew = u - h * x / rCubed(x, y)
        vNew = v - h * y / rCubed(x, y)
        xValues.append(xNew)
        yValues.append(yNew)
        x, y, u, v = xNew, yNew, uNew, vNew
    plt.xlim([-1.5, 3.3])
    plt.ylim([-1.8, 1.8])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(xValues, yValues, '-g', model(e)[0], model(e)[1], ':k', 0, 0, 'ok')
    
def FE2(d, e, N, I):
    x, y, u, v = 1 + e, 0, 0, math.sqrt((1 - e) / (1 + e))
    x0, v0 = x, v
    h = 2 * math.pi / d
    for i in range(I):
        for j in range(d * N):
            xNew, yNew = x + h * u, y + h * v
            uNew = u - h * x / rCubed(x, y)
            vNew = v - h * y / rCubed(x, y)
            x, y, u, v = xNew, yNew, uNew, vNew
        u, v = -u, -v
    return (abs(x - x0), abs(y), abs(u), abs(v - v0))  

    