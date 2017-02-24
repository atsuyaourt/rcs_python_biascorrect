import os
import json

import numpy as np
import pandas as pd
import xarray as xr

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage


def load_r_script_func(file_name):
    module_name = os.path.splitext(os.path.basename(file_name))[0]
    with open(file_name, 'r') as myfile:
        string = myfile.read()

    return SignatureTranslatedAnonymousPackage(string, module_name)


def get_aphro_varname(var_name):
    _aphro_var = {'tas': 'tave', 'pr': 'precip'}
    return _aphro_var[var_name]


def get_aphro_input(var_name, input_dir='input/aphro'):
    return ('%s/aphrodite_%s.nc'
            % (input_dir, get_aphro_varname(var_name)))


def get_aphro_output(var_name, output_dir='output/aphro'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return '%s/aphrodite_%s.nc' % (output_dir, var_name)


def get_mod_input(var_name, exp_name, gcm_name, input_dir='input/mod'):
    if var_name == 'pr':
        return '%s/%s_%s_pr.nc' % (input_dir, exp_name, gcm_name)
    return '%s/%s_%s.nc' % (input_dir, exp_name, gcm_name)


def get_mod_output(var_name, exp_name, gcm_name, output_dir='output/mod'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return '%s/%s_%s_%s.nc' % (output_dir, exp_name, gcm_name, var_name)


def get_bc_output(var_name, exp_name, gcm_name, output_dir='output/bc'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return '%s/%s_%s_%s.nc' % (output_dir, exp_name, gcm_name, var_name)


def xr_to_df(filename, var, lon_bnds, lat_bnds, time_bnds):
    ds = xr.open_dataset(filename)

    if isinstance(var, dict):
        var_name = var['name']
    elif isinstance(var, str):
        var_name = var

    dat = ds[var_name].sel(
        lat=slice(*lat_bnds),
        lon=slice(*lon_bnds),
        time=slice(*time_bnds)
    )

    if isinstance(var, dict):
        if 'add' in var:
            dat += var['add']
        elif 'mult' in var:
            dat *= var['mult']

    dat['time'].values = np.array([
        np.datetime64(dt[:10]) for dt in np.datetime_as_string(
            dat['time'].values,
            timezone='UTC')])
    return dat.to_dataframe().reset_index()


def do_qmap(df, fit_func):
    _lat = df['lat'].values[0]
    _lon = df['lon'].values[0]
    _month = df['month'].values[0]
    fit = (
        fit_func[
            (fit_func['lat'] == _lat) &
            (fit_func['lon'] == _lon) &
            (fit_func['month'] == _month)
        ].values[0][3]
    )
    if fit:
        ret_df = pandas2ri.ri2py(bias_correct.do_qmap(df, fit))
    else:
        ret_df = pandas2ri.ri2py(bias_correct.do_qmap(df))

    ret_df['time'] = df['time'].values
    return ret_df


def df_to_nc(df, out_file):
    var_name = df.name
    out_dat = df.to_xarray().to_dataset()
    vars_attr = {
        'pr': {
            'standard_name': 'precipitation',
            'long_name': 'Precipitation',
            'units': 'mm'
        },
        'tas': {
            'standard_name': 'temperature',
            'long_name': 'Temperature',
            'units': 'C'
        }
    }
    out_dat[var_name].attrs = vars_attr[var_name]
    out_dat['lat'].attrs = {
        'standard_name': 'latitude',
        'long_name': 'Latitude',
        'units': 'degrees_north'
    }
    out_dat['lon'].attrs = {
        'standard_name': 'longitude',
        'long_name': 'Longitude',
        'units': 'degrees_east'
    }
    out_dat.to_netcdf(out_file, mode='w', unlimited_dims=['time'])


pandas2ri.activate()
bias_correct = load_r_script_func('scripts/bias_correct.R')

with open('settings.json') as json_data:
    conf = json.load(json_data)

for var in conf['variables']:
    # Load APHRODITE data as pandas DataFrame
    obs_file = get_aphro_input(var['name'])
    print 'Reading %s...' % obs_file
    obs_time_bnds = ['1971-01-01', '2000-12-31']
    aphro_var = get_aphro_varname(var['name'])
    obs_df = (
        xr_to_df(obs_file, aphro_var,
                 conf['lon_bnds'], conf['lat_bnds'], obs_time_bnds)
        .drop(['lev'], axis=1)
    )

    # Save as netcdf
    out_df = obs_df.set_index(['time', 'lat', 'lon'])[aphro_var]
    out_df.name = var['name']
    out_file = get_aphro_output(var['name'])
    df_to_nc(out_df, out_file)
    print 'Writing %s done!!!' % out_file

    for gcm in conf['gcms']:
        fit_func = pd.DataFrame()
        for exp in conf['experiments']:
            # Load data as pandas DataFrame
            mod_file = get_mod_input(var['name'], exp['name'], gcm)
            print 'Reading %s...' % mod_file
            mod_df = xr_to_df(
                mod_file, var,
                conf['lon_bnds'], conf['lat_bnds'], exp['time_bnds']
            )

            if gcm == 'HadGEM2':
                mod_df = (
                    mod_df.groupby(['lon', 'lat', 'time'])
                    .mean().reset_index()
                )

            # Save as netcdf
            out_df = mod_df.set_index(['time', 'lat', 'lon'])[var['name']]
            out_file = get_mod_output(var['name'], exp['name'], gcm)
            df_to_nc(out_df, out_file)
            print 'Writing %s done!!!' % out_file

            # Fit quantile maps
            if exp['name'] == 'hist':
                print 'Fitting quantile map'
                df = obs_df.merge(
                    mod_df,
                    how='right',
                    on=['lon', 'lat', 'time']
                )
                df.rename(columns={
                    'tave': 'obs', 'tas': 'mod',
                    'precip': 'obs', 'pr': 'mod'
                }, inplace=True)
                df['month'] = pd.DatetimeIndex(df['time']).month
                fit_func = (
                    df
                    .groupby(
                        ['lon', 'lat', 'month'],
                        sort=False, as_index=False
                    )
                    .apply(bias_correct.fit_qmap)
                    .reset_index()
                )
            else:
                df = mod_df
                df.rename(columns={'tas': 'mod', 'pr': 'mod'}, inplace=True)
                df['month'] = pd.DatetimeIndex(df['time']).month

            print 'Applying bias correction...'
            out_df = (
                df.groupby(['lon', 'lat', 'month'], sort=False, as_index=False)
                .apply(do_qmap, fit_func)
                .set_index(['time', 'lat', 'lon'])
            )
            out_df = out_df['bc']
            out_df.name = var['name']
            out_file = get_bc_output(var['name'], exp['name'], gcm)
            df_to_nc(out_df, out_file)
            print 'Writing %s done!!!' % out_file
