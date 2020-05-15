import logging

logger = logging.getLogger('mstools')
logger.setLevel(logging.INFO)
# formatter = logging.Formatter('%(asctime)s    %(levelname)-10s %(message)s')
formatter = logging.Formatter('%(levelname)-10s %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)
