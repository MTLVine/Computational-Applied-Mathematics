import math
import matplotlib.pyplot as plt

def rCubed(x, y):
    return pow(pow(x, 2) + pow(y, 2), 3 / 2)

def dt(x, y, u, v):
    return (u, v, -x / rCubed(x, y), -y / rCubed(x, y))

def L(x, y, u, v):
    return x * v - y * u

def momentum(d, e, N):
    x, y, u, v = 1 + e, 0, 0, math.sqrt((1 - e) / (1 + e))
    h = 2 * math.pi / d
    L0 = L(x, y, u, v)
    time, FE, ME, LF, RK = [], [0], [0], [0], [0]
    for i in range(d * N):
        time.append(i * h)
    while len(FE) < d * N:
        xNew, yNew = x + h * u, y + h * v
        uNew = u - h * x / rCubed(x, y)
        vNew = v - h * y / rCubed(x, y)
        numerator = abs(L(xNew, yNew, uNew, vNew) - L0)
        FE.append(numerator / abs(L0))
        x, y, u, v = xNew, yNew, uNew, vNew
    x, y, u, v = 1 + e, 0, 0, math.sqrt((1 - e) / (1 + e))
    while len(LF) < d * N:
        X, Y = x + u * h / 2, y + v * h / 2
        uNew = u - h * X / rCubed(X, Y)
        vNew = v - h * Y / rCubed(X, Y)
        xNew, yNew = X + uNew * h / 2, Y + vNew * h / 2
        numerator = abs(L(xNew, yNew, uNew, vNew) - L0)
        LF.append(numerator / abs(L0))
        x, y, u, v = xNew, yNew, uNew, vNew
    x, y, u, v = 1 + e, 0, 0, math.sqrt((1 - e) / (1 + e))    
    while len(ME) < d * N:
        xNew, yNew = x + h * u, y + h * v
        uNew = u - h * xNew / rCubed(xNew, yNew)
        vNew = v - h * yNew / rCubed(xNew, yNew)
        numerator = abs(L(xNew, yNew, uNew, vNew) - L0)
        ME.append(numerator / abs(L0))
        x, y, u, v = xNew, yNew, uNew, vNew
    x, y, u, v, h = 1 + e, 0, 0, math.sqrt((1 - e) / (1 + e)), 4 * h
    while len(RK) < d * N:
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
        numerator = abs(L(x, y, u, v) - L0)
        RK.append(numerator / abs(L0))
    plt.xlabel('Time')
    plt.ylabel('Fractional Error in Angular Momentum')
    plt.loglog(time, FE, ':g', time, ME, '-g', time, LF, '-k', time, RK, ':k')