#!/usr/bin/env python

'''
This script can be used to concatenate multiple files of the same variable
for CMIP output, i.e., when there are multiple files for one experiment
split over time. A single output file is created. By default, the original
files will be deleted, but they can be kept if the -k option is used.

Author: H. Frenzel, School of Oceanography, University of Washington
(hfrenzel@uw.edu)
First version: March 2, 2021
Latest revision: January 10, 2022
'''

import argparse
import glob
import os
import pathlib
import re

def parse_filename(file, project):
    '''Parse the given name of a file, detect and return its components:
    variable, frequency, model, experiment, variant, grid, start, end
    For CMIP5 (project==5), there is no grid variable (it will be set to None).
    Returns None values if the name does not match the expected pattern.'''
    # Assumption: model names use hyphens, not underscores
    # Example file name:
    # tas_Amon_MPI-ESM1-2-LR_ssp126_r1i1p1f1_gn_201501-203412.nc
    if project == 6:
        fname_regex = re.compile(r'([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_' +
                                 r'([a-z0-9]+)_(\d+)-(\d+).nc')
        match_obj = fname_regex.search(file)
        if match_obj:
            return match_obj.groups()
    else:   # project == 5
        fname_regex = re.compile(r'([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_(\d+)-(\d+).nc')
        match_obj = fname_regex.search(file)
        if match_obj:
            variable, freq, model, experiment, variant, start, end = \
                match_obj.groups()
            return variable, freq, model, experiment, variant, None, start, end
    return None, None, None, None, None, None, None, None


def get_types(pattern, project):
    '''Create an inventory of types (combinations of variables, frequencies,
    models, experiments, variants, grids) in the form of a dictionary,
    where the values are lists with the names of the files of the given type.
    Return the dictionary. (If there are no netcdf files present in the
    current directory, the dictionary will be empty.)
    Note that the files within each list are not sorted.'''
    types = dict()
    if isinstance(pattern, str):
        all_files = glob.glob(pattern)
    else:
        # if the pattern is specified as argument, the shell expands it into
        # a list of files already
        all_files = pattern
    for file in all_files:
        variable, freq, model, experiment, variant, grid = \
            parse_filename(file, project)[0:6]
        if variable:
            if (variable, freq, model, experiment, variant, grid) not in types:
                types[variable, freq, model, experiment, variant, grid] = list()
            types[variable, freq, model, experiment, variant, grid].append(file)

    return types


def sort_files_in_types(types, project):
    '''For each type in dictionary types, sort all of its files
    by their start dates. Return the dictionary with sorted file lists.'''
    types_sorted = dict()
    for type1 in types.keys():
        all_files = types[type1]
        if len(all_files) == 1:
            types_sorted[type1] = all_files
        else:
            all_starts = list()
            for file in all_files:
                variable, start = parse_filename(file, project)[0:7:6]
                if variable:
                    all_starts.append((file, int(start)))
            sorted_files = sorted(all_starts, key=lambda st: st[1])
            new_files = list()
            for file in sorted_files:
                new_files.append(file[0])
            types_sorted[type1] = new_files
    return types_sorted


def check_file_sequence(files, project):
    '''Check if all files fit together in time. If they do (and there
    is more than one of them), return True. Return False in all other
    cases.'''
    if not files: # this should not happen!
        print('No files found')
        return False
    if len(files) == 1:
        print('Only one file found, nothing to do')
        return False

    last_end = None
    for file in files:
        variable, freq, model, exp, var, grid, start, end = \
            parse_filename(file, project)
        if not variable:
            print('Parsing error: {0:s}'.format(file))
            return False

        if last_end:
            if 'Omon' in file or 'Amon' in file:
                # format for start and end is: YYYYMM
                # for most projects, periods run from Jan to Dec,
                # but some Hadley files for CMIP5 run from Dec to Nov
                if last_end % 100 == 12:
                    offset = 89
                else:
                    offset = 1
            elif 'Oyr' in file or 'Ayr' in file:
                # format for start and end is: YYYY
                offset = 1
            else:
                print('Unknown date format')
                return False
            if int(start) != last_end + offset:
                print('Break in sequence at {0:s}'.format(file))
                return False
            last_end = int(end)
    return True


def create_ncrcat_command(files):
    '''Define and return the ncrcat command to create the concatenated file.
    Also return the name of the resulting file.'''
    variable0, freq0, model0, experiment0, variant0, grid0, start0 = \
        parse_filename(files[0], ARGS.project)[0:7]
    end1 = parse_filename(files[-1], ARGS.project)[7]
    if grid0:
        file_out = '{0:s}_{1:s}_{2:s}_'.format(variable0, freq0, model0) + \
            '{0:s}_{1:s}_{2:s}_'.format(experiment0, variant0, grid0) + \
            '{0:s}-{1:s}.nc'.format(start0, end1)
    else:
        file_out = '{0:s}_{1:s}_{2:s}_'.format(variable0, freq0, model0) + \
            '{0:s}_{1:s}_'.format(experiment0, variant0) + \
            '{0:s}-{1:s}.nc'.format(start0, end1)

    cmd = 'ncrcat '
    if ARGS.overwrite:
        cmd += '-O '
    for file in files:
        cmd += '{0:s} '.format(file)
    cmd += file_out
    return cmd, file_out


def delete_files(files, test_only):
    '''Delete all the given files without further asking, unless test_only
    is True (only print names of files that would be deleted in that case).'''
    for file in files:
        print('deleting {0:s}'.format(file))
        if not test_only:
            os.remove(file)


def check_sequence(files, project):
    '''Check if the given files form a sequence without gaps in time.
    File names depend on the given project.
    New list of files (depending on options, original files only,
    original files plus concatenated file, or concatenated file only)
    and a failure flag (for the ncrcat process) is returned.'''
    is_sequence, filename_out, failure = \
        process_file_sequence(files, ARGS.project)
    if is_sequence:
        if ARGS.keep:
            files.append(filename_out)
        else:
            files = list()
            files.append(filename_out)
    return files, failure


def process_file_sequence(files, project):
    '''Check if the list of files represents a continuous sequence in time.
    If so, concatenate the files and return True and
    the name of the resulting file.
    Otherwise return False, None.
    The third return value is True if it is a sequence, but concatenation
    failed (e.g., due to presence of zero-sized files). It is False in all
    other cases.'''
    if check_file_sequence(files, project):
        cmd, file_out = create_ncrcat_command(files)
        print(cmd)
        if not ARGS.test:
            print('concatenating, please be patient...')
            status = os.system(cmd) # 0 indicates success
            if status:
                print('concatenating failed!')
                return False, None, True
            if not ARGS.keep:
                delete_files(files, ARGS.test)
        return True, file_out, False
    return False, None, False


def move_files(this_type1, files):
    '''Move the given files to the appropriate directory. Create that
    if necessary.'''
    # directory: model/experiment/native_grid
    if not files:
        return
    dir1 = '{0:s}/{1:s}'.format(this_type1[2], this_type1[3])
    dir2 = '{0:s}/native_grid'.format(dir1)
    if not os.path.isdir(dir2):
        print('mkdir -p {0:s}'.format(dir2))
        if not ARGS.test:
            pathlib.Path(dir2).mkdir(parents=True)
            try:
                os.chmod(this_type1[2], 0o775) # drwxrwxr-x
            except OSError:
                print('Could not change permissions for {0:s}'.format(this_type1[2]))
            try:
                os.chmod(dir1, 0o775)
            except OSError:
                print('Could not change permissions for {0:s}'.format(dir1))
            try:
                os.chmod(dir2, 0o775)
            except OSError:
                print('Could not change permissions for {0:s}'.format(dir2))
    print('MOVING these files to {0:s}:'.format(dir2))
    if ARGS.test:
        print(files) # the concatenated file does not exist yet
    else:
        for file in files:
            if os.path.getsize(file):
                print(file)
                new_file = '{0:s}/{1:s}'.format(dir2, file)
                try:
                    os.rename(file, new_file)
                except OSError:
                    print('Could not move {0:s} to {1:s}'.format(file, new_file))


# https://stackoverflow.com/questions/40644096/argparse-differentiate-between-no-options-option-invoked-and-option-invoked
def parse_input_args():
    '''Parse the command line arguments and return them as object.'''
    parser = argparse.ArgumentParser(description=
                                     'Concatenate (if needed) and move CMIP ' +
                                     'files to model directories')
    # the pattern is optional - if not specified, all *nc files will be
    # processed
    parser.add_argument('pattern', nargs='*', default='*.nc')
    # option:
    parser.add_argument('-k', '--keep', default=False, action='store_true',
                        help='keep individual files')
    parser.add_argument('-n', '--no_ncrcat', default=False, action='store_true',
                        help='do not ncrcat files in sequences (move only)')
    parser.add_argument('-O', '--overwrite', default=False, action='store_true',
                        help='overwrite files if they exist already')
    parser.add_argument('-p', '--project', default=6, type=int,
                        help='project: default is CMIP6 (6); alternative is CMIP5 (5)')
    parser.add_argument('-t', '--test', default=False, action='store_true',
                        help='test only; show commands without executing them')
    args = parser.parse_args()
    if args.project != 5 and args.project != 6:
        raise ValueError('Unknown project, only 5 or 6 are allowed')
    return args


if __name__ == '__main__':
    ARGS = parse_input_args()
    all_types = get_types(ARGS.pattern, ARGS.project)
    all_types = sort_files_in_types(all_types, ARGS.project)
    for this_type in all_types:
        if not ARGS.no_ncrcat:
            all_types[this_type], failed = check_sequence(all_types[this_type],
                                                          ARGS.project)
        if ARGS.no_ncrcat or not failed:
            move_files(this_type, all_types[this_type])
