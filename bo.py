from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import copy
import sympy
import random
import time
 
class Profiler(object):
    def __enter__(self):
        self._startTime = time.time()
         
    def __exit__(self, type, value, traceback):
        print ("Elapsed time: {:.3f} sec".format(time.time() - self._startTime))


def func(args, m=10):
    # [0; pi], [-pi; pi] - область определения
    f = 0
    for i in range(len(args)):
        f += np.sin(args[i]) * np.sin((i + 1) * (args[i] ** 2)/np.pi) ** (2 * m)
    return -f

def makeData(minimum_num, maximum_num):
    # Создаем двумерную матрицу-сетку
    
    # В узлах рассчитываем значение функции
    #zgrid = -np.sin (xgrid) * np.sin (ygrid) / (xgrid * ygrid)
    zgrid = (xgrid * xgrid + ygrid * ygrid)**0.5 + 3 * np.cos((xgrid * xgrid + ygrid * ygrid)**0.5) + 5
    return xgrid, ygrid, zgrid

def generate_vector(min_length=-1, max_length=1, problem_size=3):
    size = problem_size 
    vec = np.random.uniform(min_length, max_length, size)
    return vec

def draw_3D_graphic(cells_coord, minimum_num, maximum_num):
    try:
        fig = plt.figure()
        ax = Axes3D(fig)
        x = np.arange(minimum_num, maximum_num, 0.1)
        y = np.arange(minimum_num, maximum_num, 0.1)
        xgrid, ygrid = np.meshgrid(x, y)
        zgrid = func([xgrid, ygrid])
        ax.plot_surface(xgrid, ygrid, zgrid, cmap = cm.Pastel1)
        for cell in cells_coord:
            ax.scatter(cell[0], cell[1], cell[2], s=40, c='0.1')
        ax.legend()
        #plt.show()
    except:
        print("Error")

def draw_2D_graphic(first_p, second_p, label='graph', file='figure'):
    plt.plot(first_p, second_p, label = label, linewidth = 0.8)
    plt.legend()
    # ax = plt.gca()
    # ax.set_xlim([1, 10])
    plt.savefig(file + '.png', format='png', dpi=100)
    #plt.show()
    #plt.close()


class Cell(object):
    """docstring for Cell"""
    def __init__(self, d_attr, w_attr, h_rep, w_rep):
        self.coord = None
        self.vector = None
        self.d_attr = d_attr
        self.w_attr = w_attr
        self.h_rep = h_rep
        self.w_rep = w_rep
        self.fitness = 0
        self.interaction = 0
        #self.search_space = search_space
    
    def initialize_coord(self,problem_size, min_num, max_num, func):
        """
        limit = x_min, x_max, y_min....
        """
        size = problem_size
        xy = np.random.uniform(min_num, max_num, size)
        z = func(xy)
        self.coord = np.hstack((xy, z))
        return self.coord

    def compute_cell_interaction(self, cells):
        for j in range(len(cells)):
            self.interaction += (-self.d_attr * np.exp(-self.w_attr * np.sum((self.coord - cells[j].coord) ** 2))) + (
                self.h_rep * np.exp(-self.w_rep * np.sum((self.coord - cells[j].coord) ** 2)))
        return self.interaction

    def tumble(self, step_size, vec, func):
        crd = self.coord[:-1] + step_size * vec
        self.coord = np.hstack((crd, func(crd)))


    def swim(self, step_size, func):
        vec = generate_vector()
        crd = []
        crd = self.coord[:-1] + step_size * vec
        crd = np.append(crd, func(crd[:-1]))

    def calculate_fitness(self):
        self.fitness = self.coord[-1] + self.interaction
        return self.fitness

    def __lt__(self, other):
        return self.fitness < other.fitness

    def compute_cell_interaction(self, cells):
        """ A bacteria cost is derated by its interaction with other cells. This function is calculated interaction
        """
        self.interaction = 0
        for i in range(len(cells)):
            self.interaction += (-self.d_attr * np.exp(-self.w_attr * np.sum((self.coord - cells[i].coord) ** 2)))
            + (self.h_rep * np.exp(-self.w_rep * np.sum((self.coord - cells[i].coord) ** 2)))
        return self.interaction

def chemotaxis(cells, chem_steps, swim_length, step_size,
                d_attr, w_attr, h_rep, w_rep, problem_size):
    best_cell = None
    for cs in range(chem_steps):
        for i in range(len(cells)):
            #cells[i].compute_cell_interaction(cells)
            vec = generate_vector(problem_size=problem_size)
            new_cell = copy.copy(cells[i])
            for sl in range(swim_length):
                new_cell.tumble(step_size, vec, func)
                new_cell.compute_cell_interaction(cells[:i] + cells[i+1:])
                if new_cell.calculate_fitness() <= cells[i].fitness:
                    cells[i] = copy.copy(new_cell)
                else:
                    vec = generate_vector(problem_size=problem_size)
            if best_cell == None or best_cell.fitness > cells[i].fitness:
                best_cell = copy.copy(cells[i])
            del new_cell
        #print('chemo={}, f={}, cost[{}]={}'.format(cs, best_cell.fitness, number_cell, cells[number_cell].coord))
    return best_cell, cells

def search_optimum(problem_size, cells_num, N_ed, N_re, N_c, N_s, 
    d_attract, w_attract, h_repellant, w_repellant, P_ed,
    step_size, min_num=0, max_num=5, file='results.txt', label = None):
    """ p: Dimension of the search space,
         S: Total number of bacteria in the population,
         Nc : The number of chemotactic steps,
         Ns : The swimming length.
         Nre : The number of reproduction steps,
         Ned : The number of elimination-dispersal events,
         Ped : Elimination-dispersal probability,
         C (i): The size of the step taken in the random direc
    """
    #f = open(file, 'tw', encoding='utf-8')
    r = []
    cells = []
    for i in range(cells_num):
        cells.append(Cell(d_attract, w_attract, h_repellant, w_repellant))
        cells[i].initialize_coord(problem_size, min_num, max_num, func)
    best = None
    with Profiler() as optimum_profiler:
        for ed in range(N_ed):
            print("elim_disp_steps = ", ed)
            #f.write("elim_disp_steps = {}\n".format(ed))
            for re in range(N_re):
                c_best, cells = chemotaxis(cells, N_c, N_s, step_size, d_attract, w_attract, h_repellant, w_repellant, problem_size)
                if best == None or best.fitness > c_best.fitness:
                    best = copy.copy(c_best)
                del c_best
                print('best fitness = {}, coord = {}'.format(best.fitness, best.coord))
                #f.write('best fitness = {}, coord = {} \n'.format(best.fitness, best.coord))
                cells.sort()
                cells = cells[:round(len(cells)/2)] + cells[:round(len(cells)/2)]
                r.append(best.coord[-1])
            for i in range(len(cells)):
                if random.random() < P_ed:
                    cells[i] = Cell(d_attract, w_attract, h_repellant, w_repellant)
                    cells[i].initialize_coord(problem_size, min_num, max_num, func)
                    cells[i].compute_cell_interaction(cells)
                    cells[i].calculate_fitness()
    #f.close()
    if label == None:
        label =  'Fitness for number of elimination-dispersal step'
    draw_2D_graphic(np.linspace(1, N_ed, len(r)), np.array(r), label = label, file=file[:-4])
    result = best.coord[-1]

if __name__ == '__main__':

    #default coefficients
    number_cells = 50
    d_attr = 0.1
    w_attr = 0.2
    h_repell = 0.1 # h_repell = d_attr
    w_repell = 10
    step_size = 0.01 #C[i]
    elim_disp_steps = 60 # Ned is the number of elimination-dispersal steps
    repro_steps = 4 # Nre is the number of reproduction steps
    chem_steps = 80 # Nc is the number of chemotaxis steps
    swim_length = 4 # Ns is the number of swim steps for a given cell
    p_eliminate = 0.25 # Ped

    coord = [2.0890882, 1.56956964, 1.30172081, 1.91837155, 1.71718343]
    print('function = ', func(coord))
    #The step size is commonly a small fraction of the search space, such as 0.1.
    #for i in range(cells_init.shape[0]):
    # print('cost = ', compute_cell_interaction(cells_init, 3, d_attr, w_attr, h_repell, w_repell))
    #cells_array = np.vstack((cells_init[:3], cells_init[3 + 1:]))
    #step_limit = np.array([-5, 5])
    #print(tumble(np.array([-1.0, 1.0, 0.5]), np.array([-0.9, -0.9, 1.1]), 0.9, np.array([-10.0, 10.0])))
    # ax.legend()
    # plt.show()
    search_space =  [[-10, 10], [-10, 10]]
    problem_size = 5
    min_num = -5
    max_num = 5
    #problems_size = [2,3,4,5]
    #for problem_size in problems_size:

    xyz = [2.2029, 1.5707, 1.2850, 1.9231, 1.7205]
    print('f is ', func(xyz))

    my_cell = []
    for i in range(number_cells):
    	my_cell.append(Cell(d_attr, w_attr, h_repell, w_repell))
    	my_cell[i].initialize_coord(problem_size, min_num, max_num, func)
    cells = np.array(my_cell)

    search_optimum(problem_size, number_cells, elim_disp_steps, repro_steps, chem_steps, swim_length,d_attr, w_attr, h_repell, w_repell, p_eliminate, step_size, file = 'results_.txt', label = 'function convergence for Dimension = {}'.format(problem_size))

    fig, ax = plt.subplots()
    #n, bins, patches = ax.hist(result, 50)
    ax.set_title(r'Histogram of best fitness: N = 100')
    #plt.show()

# coord = [ 2.0708629   1.57689448  1.31537946  1.91802942  1.71699419 -4.73039794]