import logging
import os
import sys

def setup_logger(work_dir):
    log_file = os.path.join(work_dir, "autocen.log")
    logger = logging.getLogger("AutoCen")
    logger.setLevel(logging.INFO)
    logger.handlers = []

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    fh = logging.FileHandler(log_file, mode='w')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    return logger
