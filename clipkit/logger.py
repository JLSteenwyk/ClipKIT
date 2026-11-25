import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)
log_file_logger = logging.getLogger("log_file")


class _SiteLogFileHandler(logging.FileHandler):
    """Buffer site records between explicit flushes."""

    def emit(self, record):
        if self.stream is None:
            return
        try:
            self.stream.write(self.format(record) + self.terminator)
        except RecursionError:
            raise
        except Exception:
            self.handleError(record)
