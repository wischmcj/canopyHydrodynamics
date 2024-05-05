import __future__
import paramiko
from src.canhydro.global_vars import log, output_dir, data_dir


def put_file(file, sftp):
    log.info(f'trying to sftp {file}')
    if './data/output/' not in file:
        start_loc = f'./{file}'
    else:
        start_loc = f'./{file}'
    end_loc = f'/home/wischmcj/Desktop/canopyHydrodynamics/{file}'
    # end_loc = f'/code/code/canopyHydrodynamics/data//{file}'
    log.info(f'sftp locs: {start_loc}, {end_loc} ')
    try:
        sftp.put(start_loc,end_loc)
        msg = f'sftp put success'
    except FileNotFoundError as e:
        msg = f'Error putting file {file}: {e}'
        log.info(msg)
    return msg


def get_file(file, sftp):
    # start_loc = f'/home/wischmcj/Desktop/canopyHydrodynamics/data/output/{file}'
    start_loc = f'/code/code/canopyHydrodynamics/data/output/{file}'
    if '/data/output/' not in file:
        end_loc = f'./data/output/{file}'
    else:
        end_loc = f'./{file}'
    log.info(f'sftp locs: {start_loc}, {end_loc} ')
    try:
        sftp.get(start_loc,end_loc)
        msg = f'sftp get success'
    except FileNotFoundError as e:
        log.info(f'Error getting file: {e}')
    return msg

def sftp(file, get=False, dest_ip = '192.168.0.94' ):
    ssh = paramiko.SSHClient() 
    log.info('def client')
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    # ssh.connect(host, username='penguaman', password='')
    ssh.connect(dest_ip, username='wischmcj',password='')
    log.info('connected')
    sftp = ssh.open_sftp()
    log.info('sftp open')
    if get: 
        msg = get_file(file, sftp)
    else: 
        msg = put_file(file, sftp)
    log.info(msg)
    return msg
    
def get(*args, **kwargs):
    return sftp(*args, **kwargs, get = True)

def put(*args, **kwargs):
    return sftp(*args, **kwargs, get = False)