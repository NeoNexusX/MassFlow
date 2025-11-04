import logging
import os
import sys
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


class ColorFormatter(logging.Formatter):
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
        original_message = super().format(record)
        levelname_fixed = f"{record.levelname + ':':<10.10}"
        # get color
        color = self.LEVEL_COLORS.get(record.levelno, Fore.RESET)
        # replace [LEVEL] with color
        colored_level = f"{color}{levelname_fixed}{Fore.RESET}"
        formatted_record = original_message.replace(f"{record.levelname}", colored_level)

        # multilines helper
        if '\n' in record.message:
            # calculate prefix length
            first_line = formatted_record.split('\n')[0]
            # (len(colored_level)-len(levelname_fixed)) = 10
            prefix_length = len(first_line) - len(record.message.split('\n')[0])-10
            # split into sever lines
            lines = record.message.split('\n')
            # add first line
            formatted_lines = [first_line]

            # add prefix
            for line in lines[1:]:
                aligned_prefix = ' ' * prefix_length
                formatted_lines.append(f"{aligned_prefix}{line}")

            # rebuild lines
            formatted_record = '\n'.join(formatted_lines)

        return formatted_record


class FileFormatter(logging.Formatter):

    def format(self, record):
        original_message = super().format(record)
        levelname_fixed = f"{record.levelname + ':':<10.10}"
        formatted_record = original_message.replace(f"{record.levelname}", levelname_fixed)

        # multilines helper
        if '\n' in record.message:
            # calculate prefix length
            first_line = formatted_record.split('\n')[0]
            prefix_length = len(first_line) - len(record.message.split('\n')[0])
            # split into sever lines
            lines = record.message.split('\n')
            # add first line
            formatted_lines = [first_line]

            # add prefix
            for line in lines[1:]:
                aligned_prefix = ' ' * prefix_length
                formatted_lines.append(f"{aligned_prefix}{line}")

            # rebuild lines
            formatted_record = '\n'.join(formatted_lines)

        return formatted_record


class LoggerManager:
    """logger manager to control the logger"""

    def __init__(self):
        self.loggers = {}
        self.log_dir = "logs"
        self.log_level = logging.INFO
        self.max_bytes = 10 * 1024 * 1024  # 10MB
        self.backup_count = 5
        self.initialized = False

    def init_app(self, log_level="info", log_dir='logs'):
        """init logger system"""
        # set logger level
        self.log_level = LOG_LEVELS.get(log_level.lower(), logging.INFO)

        # set logger dir
        if log_dir:
            self.log_dir = log_dir

        # make sure logger dir exist
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

        self.initialized = True
        return self

    def _configure_logger(self, logger, name):
        """configure a named logger with handlers"""
        # 防止重复添加handlers
        if logger.handlers:
            return logger

        # 设置logger级别
        logger.setLevel(self.log_level)

        # 防止日志传播到根logger
        logger.propagate = False

        # add console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(self.log_level)
        console_formatter = ColorFormatter(
            '%(levelname)s%(asctime)s %(lineno)d %(name)s - %(message)s',
            datefmt='%y-%m-%d %H:%M'
        )
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        # add file handler
        file_path = os.path.join(self.log_dir, f"{name}.log")
        file_handler = RotatingFileHandler(
            file_path,
            maxBytes=self.max_bytes,
            backupCount=self.backup_count
        )
        file_handler.setLevel(self.log_level)
        file_formatter = FileFormatter(
            '%(levelname)s%(asctime)s %(lineno)d %(name)s - %(message)s'
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

        return logger

    def get_logger(self, name) -> logging.Logger:
        """get or create specify logger"""
        if not self.initialized:
            self.init_app()

        if name not in self.loggers:
            logger = logging.getLogger(name)
            logger = self._configure_logger(logger, name)
            self.loggers[name] = logger

        return self.loggers[name]

    def set_level(self, level):
        """dynamic change logger level"""
        if isinstance(level, str):
            level = LOG_LEVELS.get(level.lower(), logging.INFO)

        # update logger level
        self.log_level = level

        # update all managed loggers and their handlers
        for logger in self.loggers.values():
            logger.setLevel(level)
            # update all handlers for this logger
            for handler in logger.handlers:
                handler.setLevel(level)


# 创建日志管理器单例实例
logger_manager = LoggerManager()

LOG_LEVEL = os.getenv("LOG_LEVEL", "info")

# 初始化日志管理器
logger_manager.init_app(log_level=LOG_LEVEL, log_dir=os.getenv("LOG_DIR", "logs"))


# 获取日志器的便捷函数
def get_logger(name="massflow"):
    return logger_manager.get_logger(name)
