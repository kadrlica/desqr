#!/usr/bin/env python
"""
Interface to python logging. For more info see:
https://docs.python.org/3/howto/logging.html
"""
import logging

class SpecialFormatter(logging.Formatter):
    """Custom formatter that prefixes WARNING and ERROR messages with their level."""

    FORMATS = {
        logging.WARNING: "WARNING: %(message)s",
        logging.ERROR: "ERROR: %(message)s",
    }
    DEFAULT_FORMAT = "%(message)s"

    def format(self, record):
        self._fmt = self.FORMATS.get(record.levelno, self.DEFAULT_FORMAT)
        return super().format(record)


def get_logger(name="desqr", level=logging.INFO):
    """Create and return a configured logger instance.

    Args:
        name: The logger name.
        level: The logging level (e.g. logging.DEBUG, logging.INFO).

    Returns:
        A configured logging.Logger instance.
    """
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(SpecialFormatter())
        logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
    [setattr(logger, k, v) for k, v in logging._nameToLevel.items()]
    return logger

logger = get_logger()
