import logging

logging.captureWarnings(True)
logger = logging.getLogger('mstools')
formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
screen_handler = logging.StreamHandler()
screen_handler.setLevel(logging.DEBUG)
screen_handler.setFormatter(formatter)
logger.addHandler(screen_handler)
