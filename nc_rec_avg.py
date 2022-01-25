#!/usr/bin/env python

'''Computes successive temporal averages from an input netcdf file
and writes output to a newly created file.
A typical application is to read a file with monthly values
and create a file with annual averages.
Author: H. Frenzel, School of Oceanography, University of Washington
(hfrenzel@uw.edu)

First version: April 22, 2020
Latest revision: January 17, 2022
'''

import argparse
import os
import time
import netCDF4


def copy_global_atts(nc_in, nc_out):
    '''Copies global attributes all at once via dictionary
    from nc_in to nc_out.'''
    nc_out.setncatts(nc_in.__dict__)


def def_global_atts_hist_src(nc_out, fn_in):
    '''Defines history1 and source1 (=fn_in) global attributes in nc_out.'''
    nc_out.history1 = 'Created with nc_rec_avg.py on ' + \
        time.ctime(time.time())
    nc_out.source1 = fn_in


def copy_dims(nc_in, nc_out):
    '''Copies all dimensions from nc_in to nc_out.'''
    for name, dimension in nc_in.dimensions.items():
        nc_out.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))


def copy_vars_define(nc_in, nc_out):
    '''Copies all variable definitions including attributes from nc_in
    to nc_out.'''
    for name, variable in nc_in.variables.items():
        nc_out.createVariable(name, variable.datatype, variable.dimensions)
        # copy variable attributes all at once via dictionary
        nc_out[name].setncatts(nc_in[name].__dict__)


def copy_vars_values_not_timedep(nc_in, nc_out):
    '''Copies all time-independent variables from nc_in to nc_out.
    Returns a list of the names of all time-dependent variables.'''
    vnames_tdep = list()
    for name, variable in nc_in.variables.items():
        uses_time = False
        for dim in variable.dimensions:
            if 'time' in dim:
                uses_time = True
                vnames_tdep.append(name)
                break
        if not uses_time:
            nc_out[name][:] = nc_in[name][:]
    return vnames_tdep


def get_dim(nc_in, partial_name):
    '''Returns the length of the first dimension in nc_in whose name
    includes partial_name. None is returned if no matching dimension
    was found.'''
    for name, dimension in nc_in.dimensions.items():
        if partial_name in name:
            return dimension.size
    return None


def comp_time_average(nc_in, nc_out, name, step1, nt_avg):
    '''Reads values of variable with 'name' from nc_in for 'nt_avg' time steps
    beginning with 'step1', computes the time average and writes it to
    nc_out.'''
    ndim = len(nc_in[name].dimensions)
    if ndim == 1:
        var_in = nc_in[name][step1:step1+nt_avg]
    elif 'bnd' in name or 'bound' in name:
        # special treatment for time bounds: use first lower and last upper;
        # but first make sure that this is a bounds dimension
        if len(nc_in[name].dimensions) != 2:
            raise RuntimeError('unexpected bounds dim size')
        dim2 = nc_in[name].dimensions[1]
        if len(nc_in.dimensions[dim2]) != 2:
            raise RuntimeError('unexpected bounds dim length')
        low_bnd = nc_in[name][step1,0]
        hi_bnd = nc_in[name][step1+nt_avg-1,1]
        nc_out[name][step1//nt_avg,:] = (low_bnd, hi_bnd)
        return
    else:
        var_in = nc_in[name][step1:step1+nt_avg,:]
    var_out = var_in.mean(axis=0)
    if ndim == 1:
        nc_out[name][step1//nt_avg] = var_out
    else:
        nc_out[name][step1//nt_avg,:] = var_out


def process_file(filename_in, nt_avg):
    '''This function processes the file named "filename_in".
    A new file (with "ann" replacing the "mon" part of the input file)
    is created. Each "nt_avg" steps are averaged.'''
    if not os.path.exists(filename_in):
        print('Warning: Input file "{0:s}" not found'.format(filename_in))
        return
    print('Processing {0:s}'.format(filename_in))
    filename_out = filename_in.replace('mon', 'ann')
    nc_in = netCDF4.Dataset(filename_in)
    nc_out = netCDF4.Dataset(filename_out, 'w', format='NETCDF4')

    copy_global_atts(nc_in, nc_out)
    def_global_atts_hist_src(nc_out, filename_in)
    copy_dims(nc_in, nc_out)
    copy_vars_define(nc_in, nc_out)
    vnames_tdep = copy_vars_values_not_timedep(nc_in, nc_out)

    # files may be too large to load a whole variable into memory,
    # run outer loop over time steps
    tsteps = get_dim(nc_in, 'time')

    for step in range(0, tsteps, nt_avg):
        print(step) # shows progress
        for vname in vnames_tdep:
            comp_time_average(nc_in, nc_out, vname, step, nt_avg)

    nc_in.close()
    nc_out.close()


def parse_input_args():
    '''Parse the command line arguments and return them as an object.'''
    parser = argparse.ArgumentParser(description='Computes successive ' +
                                     'averages in time')

    # required argument:
    parser.add_argument('filename_in', nargs='+',
                        help='name(s) of the input file(s)')
    # option:
    parser.add_argument('-n', '--ntimes', default = 12, type=int,
                        help='number of time steps to average over')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    ARGS = parse_input_args()
    for file in ARGS.filename_in:
        process_file(file, ARGS.ntimes)
