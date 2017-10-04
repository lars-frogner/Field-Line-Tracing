import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import subprocess
import os


def get_grid_initials(x_points=(2, 22, 5),
                      y_points=(2, 22, 5),
                      z_points=(-13, -13, 1)):

    x_initials, y_initials, z_initials = np.meshgrid(np.linspace(*x_points),
                                                     np.linspace(*y_points),
                                                     np.linspace(*z_points))
    initials = {'x_initials': x_initials.flatten(),
                'y_initials': y_initials.flatten(),
                'z_initials': z_initials.flatten()}

    return initials


def get_random_initials(x_limits=(1.5, 22.5),
                        y_limits=(1.5, 22.5),
                        z_limits=(-13, -11),
                        n_fieldlines=25):

    x_initials = np.random.uniform(low=x_limits[0], high=x_limits[1], size=n_fieldlines)
    y_initials = np.random.uniform(low=y_limits[0], high=y_limits[1], size=n_fieldlines)
    z_initials = np.random.uniform(low=z_limits[0], high=z_limits[1], size=n_fieldlines)

    initials = {'x_initials': x_initials,
                'y_initials': y_initials,
                'z_initials': z_initials}

    return initials


class fieldline_analyser:

    def __init__(self,
                 filename='fieldlines.dat',
                 use_existing_file=False,
                 x_initials=[13],
                 y_initials=[9],
                 z_initials=[-10],
                 aux_quantities=[],
                 comp_quantities=[],
                 direction='forward',
                 scheme='rkf23',
                 recompile=False,
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
                 decoupling_z=-1.3,
                 decoupling_beta=5,
                 decoupling_rate=1.0,
                 slope_correction=1.0,
                 max_fieldline_length=False,
                 max_output_points=5000,
                 min_fieldline_length=0.0,
                 merge_both_directions=True,
                 remove_top_terminated=False):

        self.u_l = 1e8

        self.all_standard_quantities = ['x', 'y', 'z', 's', 'sh']
        self.all_aux_quantities = ['e', 'r', 'r_eff', 'Pg', 'PB', 'beta']
        self.all_comp_quantities = ['shn', 'm', 'm_eff']
        self.comp_quantity_dependencies = {'shn': [], 'm': ['r'], 'm_eff': ['r_eff']}
        self.comp_quantity_functions = {'shn': self.compute_normalized_horizontal_distance,
                                        'm': self.compute_column_mass,
                                        'm_eff': self.compute_effective_column_mass}

        self.all_quantities = self.all_standard_quantities + self.all_aux_quantities + self.all_comp_quantities

        self.quantity_descriptions = {'x': 'horizontal component (x)',
                                      'y': 'horizontal component (y)',
                                      'z': 'depth',
                                      's': 'field line distance',
                                      'sh': 'horizontal field line distance',
                                      'r': 'mass density',
                                      'r_eff': 'effective mass density',
                                      'e': 'internal energy',
                                      'Pg': 'gas pressure',
                                      'PB': 'magnetic pressure',
                                      'beta': r'plasma $\beta$',
                                      'shn': 'normalized horizontal field line distance',
                                      'm': 'column mass',
                                      'm_eff': 'effective column mass'}

        self.quantity_names = {'x': r'$x$',
                               'y': r'$y$',
                               'z': r'$z$',
                               's': r'$s$',
                               'sh': r'$s_\mathrm{h}$',
                               'r': r'$\rho$',
                               'r_eff': r'$\rho_\mathrm{eff}$',
                               'e': r'$e$',
                               'Pg': r'$P_g$',
                               'PB': r'$P_B$',
                               'beta': r'$\beta$',
                               'shn': r'$s_\mathrm{h}^\mathrm{norm}$',
                               'm': r'$m$',
                               'm_eff': r'$m_\mathrm{eff}$'}

        self.quantity_units = {'x': r'$\mathrm{cm}$',
                               'y': r'$\mathrm{cm}$',
                               'z': r'$\mathrm{cm}$',
                               's': r'$\mathrm{cm}$',
                               'sh': r'$\mathrm{cm}$',
                               'r': r'$\mathrm{g}/\mathrm{cm}^3$',
                               'r_eff': r'$\mathrm{g}/\mathrm{cm}^3$',
                               'e': r'$\mathrm{erg}/\mathrm{cm}^3$',
                               'Pg': r'$\mathrm{dyn}$',
                               'PB': r'$\mathrm{dyn}$',
                               'beta': '',
                               'shn': '',
                               'm': r'$\mathrm{g}/\mathrm{cm}^2$',
                               'm_eff': r'$\mathrm{g}/\mathrm{cm}^2$'}

        self.filename = filename
        self.use_existing_file = use_existing_file
        self.x_initials = x_initials
        self.y_initials = y_initials
        self.z_initials = z_initials
        self.direction = direction
        self.scheme = scheme
        self.recompile = recompile
        self.compile_mode = compile_mode
        self.nonuniform_output = nonuniform_output
        self.constant_stepsize = constant_stepsize
        self.interp_order = interp_order
        self.interp_bias = interp_bias
        self.abs_tolerance = abs_tolerance
        self.rel_tolerance = rel_tolerance
        self.ds_out = ds_out
        self.ds_first = ds_first
        self.error_first = error_first
        self.beta = beta
        self.safety_factor = safety_factor
        self.scale_min = scale_min
        self.scale_max = scale_max
        self.decoupling_z = decoupling_z
        self.decoupling_beta = decoupling_beta
        self.decoupling_rate = decoupling_rate
        self.slope_correction = slope_correction

        self.min_fieldline_length = min_fieldline_length
        self.merge_both_directions = merge_both_directions
        self.remove_top_terminated = remove_top_terminated

        if self.direction == 'forward':
            self.direction_val = 1.0
        elif self.direction == 'backward':
            self.direction_val = -1.0
        elif self.direction == 'both':
            self.direction_val = 0.0
        else:
            raise ValueError('Invalid direction {}'.format(self.direction))

        if self.scheme not in ['rkf23', 'rkf45']:
            raise ValueError('Invalid scheme {}'.format(self.scheme))

        self.max_output_points = (round(max_fieldline_length/self.ds_out) + 1) if max_fieldline_length else max_output_points

        self.aux_quantities = [str(name)[:5] for name in aux_quantities]
        self.comp_quantities = comp_quantities
        self.var_names = ['x', 'y', 'z'] + self.aux_quantities

        self.n_fieldlines = len(self.x_initials)
        self.n_aux = len(self.aux_quantities)

        self.has_fieldlines = False

        self.compute_fieldlines()
        self.read_fieldlines()

        self.truncate_fieldlines()

        self.compute_distances(self.fieldlines)
        self.compute_distances(self.truncated)
        self.compute_distances(self.truncations)

        self.compute_comp_quantities(self.fieldlines)
        self.compute_comp_quantities(self.truncated)
        self.compute_comp_quantities(self.truncations)

    def compute_fieldlines(self):

        if self.has_fieldlines or self.use_existing_file:
            return

        print('Computing fieldlines...')

        if self.recompile:
            subprocess.check_call(['make', 'clean'])

        subprocess.check_call(['make', self.compile_mode])

        start_coordinates = [0.0]*3*self.n_fieldlines
        start_coordinates[::3] = self.x_initials
        start_coordinates[1::3] = self.y_initials
        start_coordinates[2::3] = self.z_initials

        subprocess.check_call(['./test_tracer.x',
                               '{:23.15E}'.format(self.direction_val),
                               '{:23}'.format('T' if self.nonuniform_output else 'F'),
                               '{:23}'.format('T' if self.constant_stepsize else 'F'),
                               '{:23d}'.format(self.interp_order),
                               '{:23d}'.format(self.interp_bias),
                               '{:23.15E}'.format(self.abs_tolerance),
                               '{:23.15E}'.format(self.rel_tolerance),
                               '{:23.15E}'.format(self.ds_out),
                               '{:23.15E}'.format(self.ds_first),
                               '{:23.15E}'.format(self.error_first),
                               '{:23.15E}'.format(self.beta),
                               '{:23.15E}'.format(self.safety_factor),
                               '{:23.15E}'.format(self.scale_min),
                               '{:23.15E}'.format(self.scale_max),
                               '{:23.15E}'.format(self.decoupling_z),
                               '{:23.15E}'.format(self.decoupling_beta),
                               '{:23.15E}'.format(self.decoupling_rate),
                               '{:23.15E}'.format(self.slope_correction),
                               '{:23d}'.format(self.max_output_points),
                               '{:23d}'.format(self.n_aux)] +
                              self.aux_quantities +
                              ['{:23.15E}'.format(c) for c in start_coordinates])

        os.rename('fieldlines.dat', self.filename)

        print('Done')

    def read_fieldlines(self):

        if self.has_fieldlines:
            return

        print('Reading fieldlines...', end=' ')

        with open(self.filename, 'r', encoding='utf8') as f:

            self.fieldlines = []

            dimensions = tuple(map(int, f.readline().split()))
            self.x = np.array(list(map(float, f.readline().split())))
            self.y = np.array(list(map(float, f.readline().split())))
            self.z = np.array(list(map(float, f.readline().split())))

            self.x_min = self.x[0]
            self.y_min = self.y[0]
            self.z_min = self.z[0]

            self.x_max = self.x[-1]
            self.y_max = self.y[-1]
            self.z_max = self.z[-1]

            self.x_range = self.x_max - self.x_min
            self.y_range = self.y_max - self.y_min
            self.z_range = self.z_max - self.z_min

            self.z_bottom = self.z_max

            z_top = self.z_min + 2*self.u_l*self.ds_out

            if self.direction == 'both' and self.merge_both_directions:
                factor = 0
            elif self.direction == 'both':
                factor = 2
            else:
                factor = 1

            n_included_fieldlines = 0

            if factor == 0:

                for i in range(self.n_fieldlines):

                    fieldline_length_forward = float(f.readline().strip())
                    decoupling_index = int(f.readline().strip())
                    fieldline_forward = {name: np.array(list(map(float, f.readline().split())))
                                         for name in self.var_names}

                    fieldline_length_backward = float(f.readline().strip())
                    decoupling_index = int(f.readline().strip())
                    fieldline_backward = {name: np.array(list(map(float, f.readline().split())))
                                          for name in self.var_names}

                    if (fieldline_length_forward + fieldline_length_backward >= self.min_fieldline_length):

                        self.fieldlines.append({name: np.concatenate((fieldline_backward[name][1:][::-1], fieldline_forward[name]))
                                                for name in self.var_names})

                        n_included_fieldlines += 1

            else:

                for i in range(factor*self.n_fieldlines):

                    fieldline_length = float(f.readline().strip())
                    decoupling_index = int(f.readline().strip())
                    fieldline = {name: np.array(list(map(float, f.readline().split())))
                                 for name in self.var_names}
                    fieldline['i_dec'] = decoupling_index

                    if self.remove_top_terminated and fieldline['z'][-1] < z_top:
                        continue

                    if fieldline_length >= self.min_fieldline_length:

                        self.fieldlines.append(fieldline)

                        n_included_fieldlines += 1

            self.average_elapsed_time = float(f.readline().strip())

        self.n_fieldlines = n_included_fieldlines
        self.has_fieldlines = True

        print('Done')

    def truncate_fieldlines(self):

        if not self.has_fieldlines:
            return

        print('Truncating fieldlines...', end=' ')

        self.truncated = []
        self.truncations = []

        self.n_truncated_fieldlines = 0

        for i in range(self.n_fieldlines):

            decoupling_index = self.fieldlines[i]['i_dec']

            if decoupling_index < 0:
                continue

            self.truncated.append({name: self.fieldlines[i][name][:decoupling_index]
                                              for name in self.var_names})

            self.truncations.append({name: self.fieldlines[i][name][decoupling_index:]
                                         for name in self.var_names})

            self.n_truncated_fieldlines += 1

        print('Done')

    def compute_distances(self, fieldlines):

        if not self.has_fieldlines:
            return

        print('Computing distances...', end=' ')

        x_threshold = 0.5*self.x_range
        y_threshold = 0.5*self.y_range
        z_threshold = 0.5*self.z_range

        for i in range(len(fieldlines)):

            x = fieldlines[i]['x'].copy()
            y = fieldlines[i]['y'].copy()
            z = fieldlines[i]['z'].copy()

            dx = x[1:] - x[:-1]
            dy = y[1:] - y[:-1]

            for j in range(len(dx)):
                if dx[j] > x_threshold:
                    x[j+1:] -= self.x_range
                elif dx[j] < -x_threshold:
                    x[j+1:] += self.x_range

            for j in range(len(dy)):
                if dy[j] > y_threshold:
                    y[j+1:] -= self.y_range
                elif dy[j] < -y_threshold:
                    y[j+1:] += self.y_range

            dx = x[1:] - x[:-1]
            dy = y[1:] - y[:-1]
            dz = z[1:] - z[:-1]

            ds = np.sqrt(dx**2 + dy**2 + dz**2)
            ds_h = np.sqrt(dx**2 + dy**2)

            fieldlines[i]['ds'] = ds
            fieldlines[i]['s'] = np.concatenate((np.zeros(1), np.cumsum(ds)))
            fieldlines[i]['sh'] = np.concatenate((np.zeros(1), np.cumsum(ds_h)))

        print('Done')

    def compute_comp_quantities(self, fieldlines):

        if not self.has_fieldlines:
            return

        print('Computing quantities...', end=' ')

        for quantity in self.comp_quantities:
            self.comp_quantity_functions[quantity](fieldlines)

        print('Done')

    def compute_column_mass(self, fieldlines):

        for i in range(len(fieldlines)):

            m = np.cumsum(0.5*(fieldlines[i]['r'][1:] + fieldlines[i]['r'][:-1]*fieldlines[i]['ds']))
            fieldlines[i]['m'] = np.concatenate((m[:1], m))

    def compute_effective_column_mass(self, fieldlines):

        for i in range(len(fieldlines)):

            m = np.cumsum(0.5*(fieldlines[i]['r_eff'][1:] + fieldlines[i]['r_eff'][:-1]*fieldlines[i]['ds']))
            fieldlines[i]['m_eff'] = np.concatenate((m[:1], m))

    def compute_normalized_horizontal_distance(self, fieldlines):

        for i in range(len(fieldlines)):

            sh_range = fieldlines[i]['sh'].max() - fieldlines[i]['sh'].min()
            fieldlines[i]['shn'] = (fieldlines[i]['sh'] - fieldlines[i]['sh'].min())/sh_range

    def set_equal_aspect(self, ax):

        max_range = np.array([self.x_range, self.y_range, self.z_range]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5*(self.x_max + self.x_min)
        Yb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5*(self.y_max + self.y_min)
        Zb = 0.5*max_range*np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5*(self.z_max + self.z_min)

        for xb, yb, zb in zip(Xb, Yb, Zb):
            ax.plot([xb], [yb], [zb], 'w')

    def split_boundary_crossing_fieldline(self, fieldline):

        splitted_fieldlines = []

        ds = np.sqrt((fieldline['x'][1:] - fieldline['x'][:-1])**2 +
                     (fieldline['y'][1:] - fieldline['y'][:-1])**2 +
                     (fieldline['z'][1:] - fieldline['z'][:-1])**2)

        threshold = 10*np.mean(ds)

        n = len(ds)
        i_prev = 0

        for i in range(n):
            if ds[i] > threshold:
                splitted_fieldlines.append({name: fieldline[name][i_prev:i+1]
                                            for name in fieldline})
                i_prev = i+2

        splitted_fieldlines.append({name: fieldline[name][i_prev:]
                                    for name in fieldline})

        return splitted_fieldlines

    def plot_fieldlines_3d(self, color='k', marker='-'):

        if not self.has_fieldlines:
            return

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        if color == 'auto':

            for i in range(self.n_fieldlines):

                fieldline_parts = self.split_boundary_crossing_fieldline(self.fieldlines[i])
                line, = ax.plot(fieldline_parts[0]['x'], fieldline_parts[0]['y'], fieldline_parts[0]['z'], linestyle=marker, linewidth=1)
                color = line.get_color()
                for fieldline in fieldline_parts[1:]:
                    ax.plot(fieldline['x'], fieldline['y'], fieldline['z'], color=color, linestyle=marker, linewidth=1)

        else:

            for i in range(self.n_fieldlines):
                for fieldline in self.split_boundary_crossing_fieldline(self.fieldlines[i]):
                    ax.plot(fieldline['x'], fieldline['y'], fieldline['z'], color + marker, linewidth=1)

        self.set_equal_aspect(ax)

        ax.set_xlim(self.x[0], self.x[-1])
        ax.set_ylim(self.y[0], self.y[-1])
        ax.set_zlim(self.z[0], self.z[-1])

        ax.invert_zaxis()
        ax.invert_yaxis()

        ax.set_xlabel(r'{}$\;[${}$]$'.format(self.quantity_names['x'], self.quantity_units['x']))
        ax.set_ylabel(r'{}$\;[${}$]$'.format(self.quantity_names['y'], self.quantity_units['y']))
        ax.set_zlabel(r'{}$\;[${}$]$'.format(self.quantity_names['z'], self.quantity_units['z']))

    def plot_quantities(self,
                        x_quantity='s',
                        y_quantity='z',
                        x_limits=None,
                        y_limits=None,
                        index_limits=None,
                        plot_kwargs={},
                        x_log=False,
                        y_log=False,
                        reverse_x=False,
                        reverse_y=False,
                        extra_patch=None,
                        extra_points=None,
                        extra_plot=None,
                        title=None,
                        aspect=None,
                        fieldline_set='full',
                        savename=None):

        if fieldline_set == 'full':
            fieldlines = self.fieldlines
        elif fieldline_set == 'truncated':
            fieldlines = self.truncated
        elif fieldline_set == 'truncations':
            fieldlines = self.truncations
        else:
            raise ValueError('Invalid fieldline set {}'.format(fieldline_set))

        if not 'zorder' in plot_kwargs:
            plot_kwargs['zorder'] = 0
        if not 'linewidth' in plot_kwargs:
            plot_kwargs['linewidth'] = 0.25

        if aspect is None:
            fig = plt.figure()
        else:
            fig = plt.figure(figsize=plt.figaspect(aspect))

        ax = fig.add_subplot(111)

        if index_limits is None:
            for i in range(len(fieldlines)):
                if len(fieldlines[i][x_quantity]) > 0 and len(fieldlines[i][y_quantity]) > 0:
                    ax.plot(fieldlines[i][x_quantity],
                            fieldlines[i][y_quantity],
                            **plot_kwargs)
        elif len(index_limits) == 1:
            for i in range(len(fieldlines)):
                if len(fieldlines[i][x_quantity]) > 0 and len(fieldlines[i][y_quantity]) > 0:
                    ax.plot(fieldlines[i][x_quantity][index_limits[0]],
                            fieldlines[i][y_quantity][index_limits[0]],
                            **plot_kwargs)
        elif len(index_limits) == 2:
            for i in range(len(fieldlines)):
                if len(fieldlines[i][x_quantity]) > 0 and len(fieldlines[i][y_quantity]) > 0:
                    ax.plot(fieldlines[i][x_quantity][index_limits[0]:index_limits[1]],
                            fieldlines[i][y_quantity][index_limits[0]:index_limits[1]],
                            **plot_kwargs)

        if extra_patch is not None:
            ax.add_patch(extra_patch)

        if extra_points is not None:
            ax.scatter(extra_points[0], extra_points[1], **extra_points[2])

        if extra_plot is not None:
            ax.plot(extra_plot[0], extra_plot[1], **extra_plot[2])

        if x_limits is not None:
            ax.set_xlim(*x_limits)

        if y_limits is not None:
            ax.set_ylim(*y_limits)

        if x_log:
            ax.set_xscale('log')
        if y_log:
            ax.set_yscale('log')
        if reverse_x:
            ax.invert_xaxis()
        if reverse_y:
            ax.invert_yaxis()

        ax.set_xlabel(r'{}{}'.format(self.quantity_names[x_quantity], r'' if len(self.quantity_units[x_quantity]) == 0 else r'$\;[${}$]$'.format(self.quantity_units[x_quantity])))
        ax.set_ylabel(r'{}{}'.format(self.quantity_names[y_quantity], r'' if len(self.quantity_units[y_quantity]) == 0 else r'$\;[${}$]$'.format(self.quantity_units[y_quantity])))

        if title is None:
            ax.set_title(r'{} vs. {} ({})'.format(self.quantity_descriptions[y_quantity][0].upper() + self.quantity_descriptions[y_quantity][1:], self.quantity_descriptions[x_quantity], fieldline_set))
        else:
            ax.set_title(title)

        plt.tight_layout()

        if savename:
            fig.savefig('{}.pdf'.format(savename))


def perform_fieldline_analysis(n_fieldlines=1000,
                               decoupling_z=-1.3,
                               decoupling_beta=5.0,
                               decoupling_rate=1.0,
                               use_existing_file=0):

    initials = get_random_initials(n_fieldlines=n_fieldlines)

    an_e = fieldline_analyser(filename='fieldlines_exact.dat',
                              use_existing_file=use_existing_file,
                              ds_out=1e-2,
                              aux_quantities=['beta', 'r'],
                              comp_quantities=['m'],
                              compile_mode='fast',
                              recompile=True,
                              direction='both',
                              max_fieldline_length=100,
                              min_fieldline_length=10,
                              merge_both_directions=False,
                              remove_top_terminated=True,
                              decoupling_z=decoupling_z,
                              decoupling_beta=decoupling_beta,
                              decoupling_rate=0,
                              **initials)

    an_d = fieldline_analyser(filename='fieldlines_decoupled.dat',
                              use_existing_file=use_existing_file,
                              ds_out=1e-2,
                              aux_quantities=['beta', 'r', 'r_eff'],
                              comp_quantities=['m', 'm_eff'],
                              compile_mode='fast',
                              recompile=True,
                              direction='both',
                              max_fieldline_length=100,
                              min_fieldline_length=10,
                              merge_both_directions=False,
                              remove_top_terminated=True,
                              decoupling_z=decoupling_z,
                              decoupling_beta=decoupling_beta,
                              decoupling_rate=decoupling_rate,
                              **initials)

    print('Average elapsed time (exact) = {:g} s'.format(an_e.average_elapsed_time))
    print('Average elapsed time (decoupled) = {:g} s'.format(an_d.average_elapsed_time))

    z_limits = (-1.5e9, 3e8)
    z_reduced_limits = (-2e8, z_limits[1])
    r_limits=(1e-15, 1e-5)
    beta_limits = (1e-4, 5e5)
    sh_limits = (0, 6e9)
    m_limits = (1e-7, 5e4)
    z_trunc_limits = (-2e8, 3e8)
    m_trunc_limits = (7e-4, 1e4)

    point_size = 15

    mean_truncation_m_exact = np.mean([truncation['m'][-1] for truncation in an_e.truncations])
    mean_truncation_m_decoupled = np.mean([truncation['m'][-1] for truncation in an_d.truncations])
    mean_truncation_m_eff_decoupled = np.mean([truncation['m_eff'][-1] for truncation in an_d.truncations])

    measured_slope_correction = mean_truncation_m_exact/mean_truncation_m_decoupled - 1
    slope_correction_deviation = mean_truncation_m_exact/mean_truncation_m_eff_decoupled - 1

    print('Mesured slope correction = {:g}'.format(measured_slope_correction))
    print('Slope correction deviation = {:g}'.format(slope_correction_deviation))

    # Plot beta

    an_e.plot_quantities(x_quantity='beta',
                         y_quantity='z',
                         x_limits=beta_limits,
                         y_limits=z_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_plot=([an_e.decoupling_beta,
                                      an_e.decoupling_beta,
                                      beta_limits[1]],
                                     [z_limits[1],
                                      an_e.u_l*an_e.decoupling_z,
                                      an_e.u_l*an_e.decoupling_z],
                                     {'color': 'navy',
                                      'linewidth': 1}),
                         title=r'Plasma $\beta$ along field lines',
                         fieldline_set='full',
                         savename='beta_z_full')

    # Plot density

    an_e.plot_quantities(x_quantity='r',
                         y_quantity='z',
                         x_limits=r_limits,
                         y_limits=z_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_plot=(r_limits,
                                     [an_e.u_l*an_e.decoupling_z]*2,
                                     {'color': 'navy',
                                      'linewidth': 1}),
                         title='Mass density',
                         fieldline_set='full',
                         savename='r_z_full')

    an_e.plot_quantities(x_quantity='r',
                         y_quantity='z',
                         x_limits=(9e-11, r_limits[1]),
                         y_limits=z_reduced_limits,
                         index_limits=[-1],
                         plot_kwargs={'marker': '.'},
                         x_log=True,
                         reverse_y=True,
                         extra_plot=(r_limits,
                                     [an_e.u_l*an_e.decoupling_z]*2,
                                     {'color': 'navy',
                                      'linewidth': 1}),
                         title='Mass density at termination points',
                         fieldline_set='truncated',
                         savename='r_z_truncated')

    # Plot height profiles

    an_e.plot_quantities(x_quantity='sh',
                         y_quantity='z',
                         x_limits=sh_limits,
                         y_limits=z_limits,
                         reverse_y=True,
                         extra_plot=(sh_limits,
                                     [an_e.u_l*an_e.decoupling_z]*2,
                                     {'color': 'navy',
                                      'linewidth': 1}),
                         title='Field line height profiles',
                         fieldline_set='full',
                         savename='sh_z_full_exact')

    an_d.plot_quantities(x_quantity='sh',
                         y_quantity='z',
                         x_limits=sh_limits,
                         y_limits=z_limits,
                         reverse_y=True,
                         fieldline_set='full',
                         savename='sh_z_full_decoupled')

    # Plot column masses

    an_e.plot_quantities(x_quantity='m',
                         y_quantity='z',
                         x_limits=m_limits,
                         y_limits=z_reduced_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_points=([mean_truncation_m_decoupled,
                                        mean_truncation_m_eff_decoupled,
                                        mean_truncation_m_exact],
                                       [an_e.z_bottom]*3,
                                       {'s': [point_size]*3,
                                        'c': ['crimson', 'navy', 'forestgreen']}),
                         extra_plot=(m_limits,
                                     [an_e.u_l*an_e.decoupling_z]*2,
                                     {'color': 'navy',
                                      'linewidth': 1}),
                         title='Column masses of full field lines',
                         aspect=0.5,
                         fieldline_set='full',
                         savename='m_z_full_exact')

    an_d.plot_quantities(x_quantity='m',
                         y_quantity='z',
                         x_limits=m_limits,
                         y_limits=z_reduced_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_points=([mean_truncation_m_decoupled,
                                        mean_truncation_m_eff_decoupled,
                                        mean_truncation_m_exact],
                                       [an_d.z_bottom]*3,
                                       {'s': [point_size]*3,
                                        'c': ['crimson', 'navy', 'forestgreen']}),
                         extra_plot=(m_limits,
                                     [an_e.u_l*an_e.decoupling_z]*2,
                                     {'color': 'navy',
                                      'linewidth': 1}),
                         fieldline_set='full',
                         savename='m_z_full_decoupled')

    an_d.plot_quantities(x_quantity='m_eff',
                         y_quantity='z',
                         x_limits=m_limits,
                         y_limits=z_reduced_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_points=([mean_truncation_m_decoupled,
                                        mean_truncation_m_eff_decoupled,
                                        mean_truncation_m_exact],
                                       [an_d.z_bottom]*3,
                                       {'s': [point_size]*3,
                                        'c': ['crimson', 'navy', 'forestgreen']}),
                         extra_plot=(m_limits,
                                     [an_e.u_l*an_e.decoupling_z]*2,
                                     {'color': 'navy',
                                      'linewidth': 1}),
                         title='Corrected column masses of extended field lines',
                         aspect=0.5,
                         fieldline_set='full',
                         savename='m_eff_z_full_decoupled')

    # Plot truncated parts of column masses

    an_e.plot_quantities(x_quantity='m',
                         y_quantity='z',
                         x_limits=m_trunc_limits,
                         y_limits=z_trunc_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_points=([mean_truncation_m_decoupled,
                                        mean_truncation_m_eff_decoupled,
                                        mean_truncation_m_exact],
                                       [an_e.z_bottom]*3,
                                       {'s': [point_size]*3,
                                        'c': ['crimson', 'navy', 'forestgreen']}),
                         fieldline_set='truncations',
                         savename='m_z_truncations_exact')

    an_d.plot_quantities(x_quantity='m',
                         y_quantity='z',
                         x_limits=m_trunc_limits,
                         y_limits=z_trunc_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_points=([mean_truncation_m_decoupled,
                                        mean_truncation_m_eff_decoupled,
                                        mean_truncation_m_exact],
                                       [an_d.z_bottom]*3,
                                       {'s': [point_size]*3,
                                        'c': ['crimson', 'navy', 'forestgreen']}),
                         fieldline_set='truncations',
                         savename='m_z_truncations_decoupled')

    an_d.plot_quantities(x_quantity='m_eff',
                         y_quantity='z',
                         x_limits=m_trunc_limits,
                         y_limits=z_trunc_limits,
                         x_log=True,
                         reverse_y=True,
                         extra_points=([mean_truncation_m_decoupled,
                                        mean_truncation_m_eff_decoupled,
                                        mean_truncation_m_exact],
                                       [an_d.z_bottom]*3,
                                       {'s': [point_size]*3,
                                        'c': ['crimson', 'navy', 'forestgreen']}),
                         fieldline_set='truncations',
                         savename='m_eff_z_truncations_decoupled')


if __name__ == '__main__':

    perform_fieldline_analysis(n_fieldlines=300,
                               decoupling_z=-1.3,
                               decoupling_beta=5.0,
                               decoupling_rate=1.0,
                               use_existing_file=0)
