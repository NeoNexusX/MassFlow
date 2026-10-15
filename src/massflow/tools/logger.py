import logging
import os
import sys
import warnings
from logging.handlers import RotatingFileHandler
from colorama import init, Fore, Back, Style

init(autoreset=True)  # 自动重置颜色

# 日志级别映射
LOG_LEVELS = {
    "debug": logging.DEBUG,
    "info": logging.INFO,
    "warning": logging.WARNING,
    "error": logging.ERROR,
    "critical": logging.CRITICAL
}

LOG_FORMAT = "%(levelname)s%(asctime)s %(lineno)d %(name)s - %(message)s"
CONSOLE_DATEFMT = "%y-%m-%d %H:%M"


class FileFormatter(logging.Formatter):
    """文件日志格式化器，保持原始记录名"""

    @staticmethod
    def _level_label(record) -> str:
        return f"{record.levelname + ':':<10.10}"

    def format(self, record):
        formatted = super().format(record)
        # levelname leads LOG_FORMAT: replace only that occurrence, never a
        # level name that appears inside the message text.
        formatted = formatted.replace(record.levelname, self._level_label(record), 1)

        # multilines helper: align continuation lines with the message column;
        # exception/stack text that follows the message is kept unchanged.
        message_lines = record.message.split('\n')
        if len(message_lines) > 1:
            lines = formatted.split('\n')
            aligned_prefix = ' ' * (len(lines[0]) - len(message_lines[0]))
            lines[1:len(message_lines)] = [aligned_prefix + line for line in message_lines[1:]]
            formatted = '\n'.join(lines)

        return formatted


class ColorFormatter(FileFormatter):
    """带颜色的日志格式化器（仅对控制台生效）"""

    # 定义不同日志级别的颜色
    LEVEL_COLORS = {
        logging.DEBUG: Fore.CYAN,  # 青色
        logging.INFO: Fore.GREEN,  # 绿色
        logging.WARNING: Fore.YELLOW,  # 黄色
        logging.ERROR: Fore.RED,  # 红色
        logging.CRITICAL: Fore.RED + Back.WHITE + Style.BRIGHT,  # 红底白字
    }

    def format(self, record):
        # color the label after alignment so escape codes never shift the indent
        label = self._level_label(record)
        color = self.LEVEL_COLORS.get(record.levelno, Fore.RESET)
        return super().format(record).replace(label, f"{color}{label}{Fore.RESET}", 1)


def custom_formatwarning(message, category, filename, lineno, line=None): #pylint: disable=unused-argument
    """Custom warning formatter to show a clean message with concise origin info."""
    # 仅提取文件名，去掉冗长的绝对路径
    short_filename = os.path.basename(filename)
    return f"[{short_filename}:{lineno} | {category.__name__}] {message}"

def console_stream():
    """Console log stream: stdout, or stderr when MASSFLOW_LOG_STREAM=stderr.

    An MCP stdio server (mzorbit) sets it: its stdout is the protocol channel.
    """
    if os.environ.get("MASSFLOW_LOG_STREAM", "").strip().lower() == "stderr":
        return sys.stderr
    return sys.stdout

class LoggerManager:
    """logger manager to control the logger"""

    def __init__(self):
        self.loggers = {}
        self.log_dir = "logs"
        self.log_level = logging.INFO
        self.max_bytes = 1 * 1024 * 1024  # 1MB
        self.backup_count = 5
        self.initialized = False
        self._file_handlers: dict[str, RotatingFileHandler] = {}


    def init_app(self, log_level="info", log_dir='logs'):
        """init logger system"""
        # set logger level
        self.log_level = LOG_LEVELS.get(log_level.lower(), logging.INFO)

        # set logger dir
        if log_dir:
            self.log_dir = log_dir

        # make sure logger dir exist
        os.makedirs(self.log_dir, exist_ok=True)

        # configure root logger
        self._configure_root_logger()

        self.initialized = True
        self.get_logger("massflow").info(f"Logger initialized with level {log_level} and log directory '{self.log_dir}' ")
        return self


    def _file_handler(self, file_name: str) -> RotatingFileHandler:
        """Return the rotating handler for ``file_name``, shared by all loggers writing to it.

        Two handlers rotating the same file would rename it out from under each other.
        """
        path = os.path.abspath(os.path.join(self.log_dir, file_name))
        handler = self._file_handlers.get(path)
        if handler is None:
            handler = RotatingFileHandler(path, maxBytes=self.max_bytes, backupCount=self.backup_count)
            handler.setFormatter(FileFormatter(LOG_FORMAT))
            self._file_handlers[path] = handler
        handler.setLevel(self.log_level)
        return handler


    def _add_handlers(self, logger, file_name: str) -> None:
        """Attach a colored console handler and a rotating file handler."""
        console_handler = logging.StreamHandler(console_stream())
        console_handler.setLevel(self.log_level)
        console_handler.setFormatter(ColorFormatter(LOG_FORMAT, datefmt=CONSOLE_DATEFMT))
        logger.addHandler(console_handler)
        logger.addHandler(self._file_handler(file_name))


    def _configure_root_logger(self):
        """configure root logger"""
        # get root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(self.log_level)

        # clear logger handlers (iterate over a copy: removeHandler mutates the list)
        for handler in list(root_logger.handlers):
            root_logger.removeHandler(handler)

        self._add_handlers(root_logger, "massflow.log")

        logging.captureWarnings(True)
        warnings.formatwarning = custom_formatwarning


    def _configure_logger(self, logger, name):
        """Configure a named logger with console and rotating file handlers."""
        # Avoid adding duplicate handlers
        if logger.handlers:
            return logger

        # Set base logger level
        logger.setLevel(self.log_level)

        # Prevent propagation to root logger
        logger.propagate = False

        self._add_handlers(logger, f"{name}.log")

        return logger


    def _logger_group(self, name: str) -> str:
        """Top-level logger group: 'bootstrap.config' -> 'bootstrap'."""
        return (name or "app").split(".", 1)[0]


    def get_logger(self, name="app") -> logging.Logger:
        """Get logger.

        Desired behavior:
        - child loggers (e.g. bootstrap.xxx) print to console and write to *parent file* (bootstrap.log)
        - file output keeps original record.name (bootstrap.xxx)
        """
        if not self.initialized:
            raise RuntimeError("LoggerManager is not initialized. Call LOGGER_MANAGER.init_app() first.")

        group = self._logger_group(name)

        # 1) Ensure the parent/group logger is configured ONCE (handlers live here)
        if group not in self.loggers:
            group_logger = logging.getLogger(group)
            group_logger = self._configure_logger(group_logger, group)  # writes to logs/{group}.log
            self.loggers[group] = group_logger

        # 2) If requesting the group itself, return it
        if name == group:
            return self.loggers[group]

        # 3) Child logger: no handlers; propagate to parent/group so it uses parent's handlers
        if name not in self.loggers:
            child = logging.getLogger(name)

            # Let records bubble up to 'bootstrap' (parent), which has the file+console handlers
            child.propagate = True

            # Keep it permissive; filtering/levels are enforced by parent/group effective level
            child.setLevel(logging.NOTSET)

            self.loggers[name] = child

        return self.loggers[name]


    def set_level(self, level):
        """dynamic change logger level"""
        if isinstance(level, str):
            level = LOG_LEVELS.get(level.lower(), logging.INFO)

        # update logger level
        self.log_level = level

        # update root logger level
        root_logger = logging.getLogger()
        root_logger.setLevel(level)

        # update root handler
        for handler in root_logger.handlers:
            handler.setLevel(level)

        # update all logger
        for logger in self.loggers.values():
            logger.setLevel(level)
            # update all handlers for this logger
            for handler in logger.handlers:
                handler.setLevel(level)

# 创建日志管理器单例实例
LOGGER_MANAGER = LoggerManager()

LOGGER_MANAGER.init_app()

# 获取日志器的便捷函数
def get_logger(name="app") -> logging.Logger:
    """Get logger by name."""
    return LOGGER_MANAGER.get_logger(name)
