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

def dt(x, y, u, v):
    return (u, v, -x / rCubed(x, y), -y / rCubed(x, y))

def RK(d, e, N):
    x, y, u, v = 1 + e, 0, 0,  math.sqrt((1 - e) / (1 + e))
    h = 2 * math.pi / (d / 4)
    xValues, yValues = [x], [y]
    while len(xValues) < d * N:
        w = [x, y, u, v]
        k1 = dt(x, y, u, v)
        k1 = [h * k1[i] for i in range(4)]
        k2 = dt(x + k1[0] / 2, y + k1[1] / 2, u + k1[2] / 2, v + k1[3] / 2)
        k2 = [h * k2[i] for i in range(4)]
        k3 = dt(x + k2[0] / 2, y + k2[1] / 2, u + k2[2] / 2, v + k2[3] / 2)
        k3 = [h * k3[i] for i in range(4)]
        k4 = dt(x + k3[0], y + k3[1], u + k3[2], v + k3[3])
        k4 = [h * k4[i] for i in range(4)]
        w = [w[i] + ((k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6) for i in range(4)]
        x, y, u, v = w[0], w[1], w[2], w[3]
        xValues.append(x)
        yValues.append(y)
    #plt.xlim([-1.5, 3.3])
    #plt.ylim([-1.8, 1.8])
    plt.xlim([-0.1,0.7])
    plt.ylim([0, 0.6])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(xValues, yValues, '-g', model(e)[0], model(e)[1], ':k', 0, 0, 'ok')
    
def RK2(d, e, N, I):
    x, y, u, v = (1 + e), 0, 0,  math.sqrt((1 - e) / (1 + e))
    x0, v0 = x, v
    h = 2 * math.pi / (d / 4)
    for j in range(I):
        for j in range(d * N):
            w = [x, y, u, v]
            k1 = dt(x,y,u,v)
            k1 = [h * k1[i] for i in range(4)]
            k2 = dt(x + k1[0] / 2, y + k1[1] / 2, u + k1[2] / 2, v + k1[3] / 2)
            k2 = [h * k2[i] for i in range(4)]
            k3 = dt(x + k2[0] / 2, y + k2[1] / 2, u + k2[2] / 2, v + k2[3] / 2)
            k3 = [h * k3[i] for i in range(4)]
            k4 = dt(x + k3[0], y + k3[1], u + k3[2], v + k3[3])
            k4 = [h * k4[i] for i in range(4)]
            w = [w[i] + ((k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6) for i in range(4)]
            x, y, u, v = w[0], w[1], w[2], w[3]
        u, v = -u, -v
    return (abs(x - x0), abs(y), abs(u), abs(v - v0))