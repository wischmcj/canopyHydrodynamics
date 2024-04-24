import __future__
import paramiko
from src.canhydro.global_vars import log, output_dir, data_dir


def put_file(file, sftp):
    start_loc = f'./{file}'
    end_loc = f'/home/wischmcj/Desktop/canopyHydrodynamics/{file}'
    log.info(f'sftp locs: {start_loc}, {end_loc} ')
    try:
        put_result =sftp.put(start_loc,end_loc)
    except FileNotFoundError as e:
        log.info('Error getting file: {e}')
    msg = f'sftp put result:{put_result}'
    return msg


def get_file(file, sftp):
    start_loc = f'/home/wischmcj/Desktop/canopyHydrodynamics/data/output/{file}'
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
    paramiko.Agent()
    ssh.connect(host, username='wischmcj',password='')
    log.info('connected')
    sftp = ssh.open_sftp()
    log.info('sftp open')
    if not get: #*******REVERSED FOR TEST RUN ON PI
        msg = get_file(file, sftp)
    else: 
        msg = put_file(file, sftp)
    log.info(msg)
    return msg
    