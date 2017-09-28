import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import subprocess

last_scheme = ''


def set_equal_aspect(x, y, z, ax):

    x_min = x.min()
    y_min = y.min()
    z_min = z.min()

    x_max = x.max()
    y_max = y.max()
    z_max = z.max()

    max_range = np.array([x_max - x_min, y_max - y_min, z_max - z_min]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5*(x_max + x_min)
    Yb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5*(y_max + y_min)
    Zb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5*(z_max + z_min)

    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')


def compute_fieldlines(x_initials=[13],
                       y_initials=[9],
                       z_initials=[-10],
                       direction='forward',
                       scheme='rkf45',
                       compile_mode='debug',
                       nonuniform_output=False,
                       constant_stepsize=False,
                       interp_order=3,
                       interp_bias=0,
                       abs_tolerance=1e-6,
                       rel_tolerance=1e-6,
                       ds_out=1e-3,
                       ds_first=1e-6,
                       error_first=1e-4,
                       beta=0.0,
                       safety_factor=0.9,
                       scale_min=0.2,
                       scale_max=10.0,
                       max_fieldline_length=False,
                       max_output_points=5000,
                       aux_names=[]):

    global last_scheme

    if direction == 'forward':
        direction_val = 1.0
    elif direction == 'backward':
        direction_val = -1.0
    elif direction == 'both':
        direction_val = 0.0
    else:
        raise ValueError('Invalid direction {}'.format(direction))

    n_fieldlines = len(x_initials)
    n_aux = len(aux_names)
    aux_names = [str(name)[:5] for name in aux_names]
    var_names = ['x', 'y', 'z'] + aux_names

    if max_fieldline_length:
        max_output_points = round(max_fieldline_length/ds_out) + 1

    if scheme not in ['rkf23', 'rkf45']:
        raise ValueError('Invalid scheme {}'.format(scheme))

    if scheme != last_scheme:

        subprocess.check_call(['make', 'clean'])

        subprocess.check_call(['rm', '-f', 'makefile', '*.x'])

        subprocess.check_call(['makemake.py',
                               'test_tracer.f90',
                               'tracer.f90',
                               'tracer_base.f90',
                               'tracer_params.f90',
                               'mesh.f90',
                               'poly_interpolation.f90',
                               'stepper_{}.f90'.format(scheme.upper())])

    subprocess.check_call(['make', compile_mode])

    start_coordinates = [0.0]*3*n_fieldlines
    start_coordinates[::3] = x_initials
    start_coordinates[1::3] = y_initials
    start_coordinates[2::3] = z_initials

    subprocess.check_call(['./test_tracer.x',
                           '{:23.15E}'.format(direction_val),
                           '{:23}'.format('T' if nonuniform_output else 'F'),
                           '{:23}'.format('T' if constant_stepsize else 'F'),
                           '{:23d}'.format(interp_order),
                           '{:23d}'.format(interp_bias),
                           '{:23.15E}'.format(abs_tolerance),
                           '{:23.15E}'.format(rel_tolerance),
                           '{:23.15E}'.format(ds_out),
                           '{:23.15E}'.format(ds_first),
                           '{:23.15E}'.format(error_first),
                           '{:23.15E}'.format(beta),
                           '{:23.15E}'.format(safety_factor),
                           '{:23.15E}'.format(scale_min),
                           '{:23.15E}'.format(scale_max),
                           '{:23d}'.format(max_output_points),
                           '{:23d}'.format(n_aux)] +
                          aux_names +
                          ['{:23.15E}'.format(c) for c in start_coordinates])

    last_scheme = scheme

    return n_fieldlines, var_names, direction


def read_fieldlines(use_existing=False, min_fieldline_length=0.0, **kwargs):

    if not use_existing:
        n_fieldlines, var_names, direction = compute_fieldlines(**kwargs)

    with open('fieldlines.dat', 'r') as f:

        fieldlines = []

        dimensions = tuple(map(int, f.readline().split()))
        x = np.array(list(map(float, f.readline().split())))
        y = np.array(list(map(float, f.readline().split())))
        z = np.array(list(map(float, f.readline().split())))

        if direction == 'both':

            for i in range(n_fieldlines):

                fieldline_length_forward = float(f.readline().strip())
                fieldline_forward = {name: np.array(list(map(float, f.readline().split())))
                                     for name in var_names}

                fieldline_length_backward = float(f.readline().strip())
                fieldline_backward = {name: np.array(list(map(float, f.readline().split())))
                                      for name in var_names}

                if (fieldline_length_forward + fieldline_length_backward >= min_fieldline_length):
                    fieldlines.append(fieldline_forward)
                    fieldlines.append(fieldline_backward)

        else:

            for i in range(n_fieldlines):

                fieldline_length = float(f.readline().strip())
                fieldline = {name: np.array(list(map(float, f.readline().split())))
                             for name in var_names}

                if fieldline_length >= min_fieldline_length:
                    fieldlines.append(fieldline)

    xs, xe, ys, ye, zs, ze, n_output_points = dimensions

    mx = xe - xs + 1
    my = ye - ys + 1
    mz = ze - zs + 1

    return x, y, z, fieldlines


def split_fieldline(fieldline):

    splitted_fieldlines = []

    dl = np.sqrt((fieldline['x'][1:] - fieldline['x'][:-1])**2 +
                 (fieldline['y'][1:] - fieldline['y'][:-1])**2 +
                 (fieldline['z'][1:] - fieldline['z'][:-1])**2)

    threshold = 10*np.mean(dl)

    n = len(dl)
    i_prev = 0

    for i in range(n):
        if dl[i] > threshold:
            splitted_fieldlines.append({name: fieldline[name][i_prev:i+1]
                                        for name in fieldline})
            i_prev = i+2

    splitted_fieldlines.append({name: fieldline[name][i_prev:]
                                for name in fieldline})

    return splitted_fieldlines


def plot_fieldlines(color='k', marker='-', **kwargs):

    x, y, z, fieldlines = read_fieldlines(**kwargs)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i in range(len(fieldlines)):
        for fieldline in split_fieldline(fieldlines[i]):
            ax.plot(fieldline['x'], fieldline['y'], fieldline['z'], color + marker, linewidth=1)

    set_equal_aspect(x, y, z, ax)

    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])
    ax.set_zlim(z[0], z[-1])

    ax.invert_zaxis()
    ax.invert_yaxis()


def plot_fieldlines_from_grid(x_points=(2, 22, 5),
                              y_points=(2, 22, 5),
                              z_points=(-13, -6, 1), **kwargs):

    x_initials, y_initials, z_initials = np.meshgrid(np.linspace(*x_points),
                                                     np.linspace(*y_points),
                                                     np.linspace(*z_points))
    x_initials = x_initials.flatten()
    y_initials = y_initials.flatten()
    z_initials = z_initials.flatten()

    plot_fieldlines(x_initials=x_initials,
                    y_initials=y_initials,
                    z_initials=z_initials, **kwargs)


if __name__ == '__main__':

    plot_fieldlines_from_grid(use_existing=0,
                              direction='both',
                              max_fieldline_length=50,
                              min_fieldline_length=15,
                              scheme='rkf45',
                              compile_mode='debug')

    plt.show()
