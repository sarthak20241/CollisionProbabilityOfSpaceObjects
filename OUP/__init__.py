import os
import logging

data_main_dir = os.path.join(os.getcwd(), "data")
date_format = "%Y-%m-%d"
date_time_format = "%Y-%m-%d %H:%M:%S'"

def change_data_path(path):
    global data_main_dir
    if os.path.isabs(path):
        data_main_dir = path
    else:
        data_main_dir = os.path.join(os.getcwd(), path)

def use_data_path(path):
    # TODO: fix to create path
    dir_path = os.path.join(data_main_dir, os.path.dirname(path))
    
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    
    path = os.path.join(data_main_dir, path)
    return path


# Create logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('oup_logger')