#! /usr/bin/env python

"""
Create and edit the source table.

"""

import numpy as np
import argparse

from astropy.table import Table, Column


def create_source_table(colname, colunit, coltype, savename, overwrite):
    """
    Creates an empty (len=1) fits table containing the columns parsed through argparse.

    parameters
    ----------
    colname : list
        List containing the column names

    colunit : list
        List containing the column units

    coltype: list
        List containing the column dtypes

    savename: str
        Name of output fits table.

    overwrite: Bool
        Overwriting permissions.
    """
    root = '/mnt/e/Thesis/Paper_II/Pipeline/'
    t = Table()
    for i in range(0,len(colname)):
        if colunit[i] == '-':
            colunit[i] = ''
        t[colname[i]] = Column(np.zeros([1]), unit=colunit[i], dtype=coltype[i])
    savepath = '{0}/Tables/{1}.fits'.format(root, savename)
    t.write(savepath, overwrite=overwrite)
    return

def create_cmp_table(colname, colunit, coltype, savename, overwrite):
    root = '/mnt/e/Thesis/Paper_II/Pipeline/'
    t = Table()
    for i in range(0,len(colname)):
        if colunit[i] == '-':
            colunit[i] = ''
        t[colname[i]] = Column(np.zeros([1]), unit=colunit[i], dtype=coltype[i])
    savepath = '{0}/Tables/{1}.fits'.format(root, savename)
    t.write(savepath, overwrite=overwrite)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prefix_chars='-')

    group1 = parser.add_argument_group("Table operation mode")
    group1.add_argument("--mode", dest='mode', type=str, default=None, help='Table operation mode. "create" to create a new table, "edit" to edit existing row, "add" to add new row.')

    group7 = parser.add_argument_group("Table type")
    group7.add_argument("--type", dest='type', type=str, default=None,help='')

    # Create table mode.
    group2 = parser.add_argument_group("Column names")
    group2.add_argument("--colname", dest='colname', nargs='+', default=None,
                        help='List containing column names.')
    group3 = parser.add_argument_group("Column units")
    group3.add_argument("--colunit", dest='colunit', nargs='+', default=None,
                        help='List containing column units.')
    group4 = parser.add_argument_group("Column data type")
    group4.add_argument("--coltype", dest='coltype', nargs='+', default=None,
                        help='List containing column data types.')
    group4.add_argument("--txt_file", dest='txt_file', type=str, default=None)
    group5 = parser.add_argument_group("Table savename")
    group5.add_argument("--savename", dest='savename', type=str, default=None,
                        help='Name of output fits table.')
    group6 = parser.add_argument_group("Overwrite table")
    group6.add_argument("--overwrite", dest='overwrite', action='store_true', default=False,
                        help='Table overwriting permissions.')

    options = parser.parse_args()

    if options.mode == 'create':
        # if options.type == 'source':
        #     create_source_table(options.colname, options.colunit, options.coltype, options.savename, options.overwrite)
        # elif options.type == 'output_cmp':
        #     create_cmp_table(options.colname, options.colunit, options.coltype, options.savename, options.overwrite)
        if options.txt_file:
            path_to_txt_file = '/mnt/e/Thesis/Paper_II/Pipeline/{0}'.format(options.txt_file)  # FIX
            print("Reading in columns from {0}".format(path_to_txt_file))
            colname_l, colunit_l, coltype_l = [], [], []
            with open(path_to_txt_file) as f:
                lines = [line.rstrip() for line in f]
                for i in range(0, len(lines)):
                    print('+ + + + + + + + + + + + + + + +')
                    colname, colunit, coltype = lines[i].split(',')
                    colname_l.append(colname.replace(' ', ''))
                    colunit_l.append(colunit.replace(' ', ''))
                    coltype_l.append(coltype.replace(' ', ''))
            print(colname_l, colunit_l, coltype_l)
            create_source_table(colname_l, colunit_l, coltype_l, options.savename, options.overwrite)

    elif options.type == 'edit':
        if options.type == 'add_row':
            pass
        elif options.type == 'add_col':
            pass
        elif options.type == 'edit_entry':
            pass
    elif options.mode == 'xmatch':
        pass