import logging

logging.captureWarnings(True)
logger = logging.getLogger('mstools')
# formatter = logging.Formatter('%(asctime)s    %(levelname)-10s %(message)s')
formatter = logging.Formatter('%(levelname)-10s %(message)s')
screen_handler = logging.StreamHandler()
screen_handler.setLevel(logging.INFO)
screen_handler.setFormatter(formatter)
logger.addHandler(screen_handler)
