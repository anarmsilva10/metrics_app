import logging
import os
from datetime import datetime

logger = None

class Logging:
    """
    Logging handler for unified logging across multiple modules.

    Notes:
    - A log file is only created if a `foldername` is stated.
    - If `foldername` is None, `console` must be True.
    - Automatically adds module name and timestamp to each log line.
    
    Args:
        basename (str, optional): Base name for the log file. Defaults to 'log'.
        foldername (str, optional): Folder to store logs. Defaults to None.
        console (bool, optional): Whether to log to console. Defaults to True.
        log_level (str, optional): Logging level (e.g., 'DEBUG', 'INFO', 'WARNING'). Defaults to 'INFO'.
    """

    def __init__(self, basename='log', foldername=None, console=True, log_level='INFO'):
        global logger
        basename, _ = os.path.splitext(basename)
        self.basename = basename
        self.foldername = foldername
        self.console = console
        self.log_level = log_level.upper()

        if self.foldername:
            os.makedirs(self.foldername, exist_ok=True)

        logger_key = "console_only" if self.foldername is None else f"{self.basename}|{self.foldername}"

        logger = logging.getLogger(logger_key)
        
        logger.setLevel(getattr(logging, self.log_level, logging.INFO))
        logger.propagate = False

        formatter = logging.Formatter(
            "%(asctime)s — %(levelname)s — %(module)s — %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        )
        if not self.console and not self.foldername:
            raise Exception("foldername can't be None if console output is set to False")
        
        if self.console:
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)

        if self.foldername:
            log_basename = f'{self.basename}_{datetime.now().strftime("%Y-%m-%d")}.log'
            log_path = os.path.join(self.foldername, log_basename)
            file_handler = logging.FileHandler(log_path)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

    def get_logger(self):
        """Returns the configured logger instance."""
        return logger

if logger is None:
	logger = Logging().get_logger()