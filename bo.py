from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def func(x, y):
    return -np.sin (x) * np.sin (y) / (x * y)

def initialize_population(cells_num, problem_size, min_num, max_num):
    size = (cells_num, problem_size)
    xy = np.random.uniform(min_num, max_num, size)
    z = func(xy[:,0], xy[:,1]).reshape((cells_num,1))
    return np.hstack((xy, z))

def BO(problem_size, cells_num, N_ed, N_re, N_e, N_s, step_size, w_attract, h_repellant, w_repellant, P_ed):
	pass

def search_optimum():
    pass

def makeData(minimum_num, maximum_num):
    x = np.arange(minimum_num, maximum_num, 0.1)
    y = np.arange(minimum_num, maximum_num, 0.1)
    # Создаем двумерную матрицу-сетку
    xgrid, ygrid = np.meshgrid(x, y)
    # В узлах рассчитываем значение функции
    zgrid = -np.sin (xgrid) * np.sin (ygrid) / (xgrid * ygrid)
    return xgrid, ygrid, zgrid



if __name__ == '__main__':
    fig = plt.figure()
    ax = Axes3D(fig)
    x, y, z = makeData(-10, 10)
    ax.plot_surface(x, y, z, color = '0.9')
    side_length = 10
    cells_init = initialize_population(20, 2, -10, 10)
    for cell in cells_init:
        ax.scatter(cell[0], cell[1], cell[2], c='0.1')
    ax.legend()
    plt.show()