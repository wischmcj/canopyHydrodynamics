import __future__
import paramiko
from src.canhydro.global_vars import log, output_dir, data_dir
import glob

def put_file(file, sftp):
    file.replace('./','')
    start_loc = f'./{file}'
    end_loc = f'/code/code/canopyHydrodynamics/{file}'
    log.info(f'{start_loc} {end_loc}')
    try:
        put_result =sftp.put(start_loc,end_loc)
        msg = f'sftp put success from {start_loc} to {end_loc}'
    except FileNotFoundError as e:
        msg =f'Error putting to {end_loc} from {start_loc} to: {e}'
    return msg


def get_file(file, sftp):
    start_loc = f'/home/wischmcj/Desktop/canopyHydrodynamics/data/output/{file}'
    end_loc = f'./{file}'
    try:
        get_result = sftp.get(start_loc,end_loc)
    except FileNotFoundError as e:
        log.info('Error getting file: {e}')
    msg = f'sftp get success from {start_loc} to {end_loc} '
    return msg

def sftp(file, get=False, host = '192.168.0.105' ):
    ssh = paramiko.SSHClient() 
    log.info(f'ssh client to {host} created')
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    # ssh.connect(host, username='penguaman', password='')
    paramiko.Agent()
    ssh.connect(host, username='penguaman',password='')
    log.info('connected')
    sftp = ssh.open_sftp()
    log.info('sftp open')
    if file:
    	to_send = [file]
    else:
        from itertools import chain
        stats_files = glob.glob('./data/output/statistics/*.csv')
        flow_files = glob.glob('./data/output/flows/*.csv')
        to_send = [x for x in chain(stats_files,flow_files)]
    log.info(f'sending files {to_send}')
    log.info(f'sending files {to_send}')
    for file in to_send:
        try:
            if get:
                msg = get_file(file, sftp)
            else:
                msg = put_file(file, sftp)
        except Exception as e:
            log.info(f'Error sftping {file}, {e}, {msg}')
    log.info(msg)

