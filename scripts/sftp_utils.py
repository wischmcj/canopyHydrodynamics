import __future__
import paramiko
from src.canhydro.global_vars import log, output_dir, data_dir


def put_file(file, sftp):
    if '/data/output/' not in file:
        start_loc = f'./data/output/{file}'
    else:
        start_loc = f'./{file}'
    # end_loc = f'/home/wischmcj/Desktop/canopyHydrodynamics/data/output/{file}'
    end_loc = f'/code/code/canopyHydrodynamics/data/output/{file}'
    log.info(f'sftp locs: {start_loc}, {end_loc} ')
    try:
        put_result =sftp.put(start_loc,end_loc)
        msg = f'sftp put result:{put_result}'
    except FileNotFoundError as e:
        msg = f'Error getting file: {e}'
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
        get_result = sftp.get(start_loc,end_loc)
    except FileNotFoundError as e:
        log.info('Error getting file: {e}')
    msg = f'sftp get result:{get_result}'
    return msg

def sftp(file, get=False, host = '192.168.0.105' ):
    ssh = paramiko.SSHClient() 
    log.info('def client')
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    # ssh.connect(host, username='penguaman', password='')
    ssh.connect(host, username='wischmcj',password='')
    log.info('connected')
    sftp = ssh.open_sftp()
    log.info('sftp open')
    if get: 
        msg = get_file(file, sftp)
    else: 
        msg = put_file(file, sftp)
    log.info(msg)
    return msg
    