from __future__ import annotations

import os
import io
import sys
import time
# import debugpy
# from random import randint

# from timeit import timeit
# from time import time 
# debugpy.listen(("0.0.0.0", 5678))
# import geopandas as geo
# import matplotlib.pyplot as plt

# import numpy as np
from flask import Flask, render_template
# import pytest

# from src.canhydro.Forester import Forester
# from src.canhydro.global_vars import log, test_input_dir
from scripts.benchmark_comp import compare, initialize_forester
from scripts.benchmark_comp import compare, initialize_forester
from test.sensitivity_analysis import sensitivity_analysis
from src.canhydro.global_vars import log, output_dir, data_dir
from scripts.sftp_utils import sftp
# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
# from matplotlib.figure import Figure
import os
import paramiko



sys.path.insert(0, os.path.dirname(os.getcwd()))

app = Flask(__name__)

@app.route('/time/<string:file>')
def time(file):
    string = compare(file)
    return  render_template('index.html', strings=[string])


@app.route('/copyOver/<string:file>')
def copyOver(file):
    ssh = paramiko.SSHClient() 
    print('def client')
    ssh.load_host_keys(os.path.expanduser(os.path.join("~", ".ssh", "known_hosts")))
    print('load_keys')
    ssh.connect('192.168.0.105', username='penguaman', password='Gamma@13')
    
    print('connect')
    sftp = ssh.open_sftp()
    print('open')
    try:
        sftp.put(f'./data/{file}', '/code/code/data_from_server/')
    except Exception as e:
        print(f'Unable to put {file}: {e}' )
    print('put')
    sftp.close()
    print('sftp closed')
    ssh.close()
    print('ssh closed')
    return  'sftp successful'


@app.route('/cases')
def run_cases():
    log.info('Starting Sensitivity Analysis...')
    string = sensitivity_analysis()

    return string

@app.route('/get/<string:file>')
def sftp_get(file):
    msg = sftp(file, get=True)
    return  msg

@app.route('/put/<string:file>')
def sftp_put(file):
    msg = sftp(file, get=False)
    return  msg



@app.route('/hello')
def hello():
    return 'Hello, World!'

# @app.route('/timesum')
# def local_run():
#     string = compare()
#     return  render_template('index.html', strings=[string])


# removed to allow removal of plotting functions 
# @app.route('/plot/<string:file>')
# def plot_png(file):
#     forest = initialize_forester(test_input_dir, file)
#     tree = forest.cylinder_collections[0]
#     tree.initialize_digraph_from()
#     tree.find_flow_components()
#     # return str(len(tree.cylinders))
#     # fig = tree.draw('YZ', highlight_lambda=lambda: is_stem, filter_lambda = lambda: cyl_id>100)
#     fig = tree.draw('YZ', highlight_lambda=lambda: is_stem, filter_lambda = lambda: cyl_id>80)
#     output = io.BytesIO()
#     FigureCanvas(fig).print_png(output)
#     return Response(output.getvalue(), mimetype='image/png')


# @app.route('/test')
# def test():
#     args_str = "test/test_collection_integration.py"
#     args = args_str.split(" ")
#     pytest.main(args)
#     return 'Tested'



if __name__ == '__main__':
    # pre_populate_cache()
    app.run(debug=True, host='0.0.0.0')


# if __name__ == "__main__":
#     forest = initialize_forester(test_input_dir,"5_SmallTree.csv")
#     tree = forest.cylinder_collections[0]
#     print("Waiting for client to attach...4")
#     min_number = int(input('Please enter the min number: '))
#     max_number = int(input('Please enter the max number: '))
#     if (max_number < min_number):
#             print('Invalid input - shutting down...')
#     else:
#         print('Thanks.')
#     debugpy.wait_for_client()
#     time.sleep(5)