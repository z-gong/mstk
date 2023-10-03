import os
import logging
from pathlib import Path

DIR_MSTK = Path(__file__).parent.as_posix()

MSTK_FORCEFIELD_PATH = [os.path.join(DIR_MSTK, 'data', 'forcefield')]
_ff_path_str = os.getenv('MSTK_FORCEFIELD_PATH')
if _ff_path_str:
    MSTK_FORCEFIELD_PATH = _ff_path_str.strip(':').split(':') + MSTK_FORCEFIELD_PATH


class _CustomFormatter(logging.Formatter):
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s %(levelname)s (%(filename)s:%(lineno)d) %(message)s"

    FORMATS = {
        logging.DEBUG   : grey + format + reset,
        logging.INFO    : grey + format + reset,
        logging.WARNING : yellow + format + reset,
        logging.ERROR   : red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%m-%d %H:%M:%S")
        return formatter.format(record)


def _get_logger():
    logger = logging.getLogger('mstk')
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(_CustomFormatter())
    logger.addHandler(handler)
    return logger


logger = _get_logger()
