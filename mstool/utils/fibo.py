import numpy as np

def fibo(r=None, N=None, APL=None, plot='plot.pdf'):
    if r is not None and N is not None:
        APL = 4 * np.pi * r**2 / N
    elif r is not None and APL is not None:
        N = 4 * np.pi * r**2 / APL
    elif N is not None and APL is not None:
        r = np.sqrt(APL * N / 4 / np.pi)
    else:
        raise ValueError('sepcify two out of three args: r, N, APL')

    N = int(N)
    print('r: %8.3f\nN: %8d\nAPL: %6.3f\n\n' %(r, N, APL))

    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(N):
        y = r * (1 - (i / (N - 1)) * 2)  # y goes from r to -r
        radius = np.sqrt(r ** 2 - y ** 2)  # radius at y

        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))

    points = np.array(points)
    if plot:

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        # aspect ratio is 1:1:1 in data space
        ax.set_box_aspect((np.ptp(points[:,0]), np.ptp(points[:,1]), np.ptp(points[:,2])))
        ax.scatter(points[:,0], points[:,1], points[:,2])
        fig.savefig(plot)

    return points

