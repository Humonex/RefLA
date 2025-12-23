#!/usr/bin/env python
# coding=utf-8
# logging模块打印log日志到文件和控制台
# Log级别：CRITICAL(50)、ERROR(40)、WARNING(30)、INFO(20)、DEBUG(10)、NOTSET(0)
# 四大组件：日志器(Logger)、处理器(Handler)、格式器(Formatter)、过滤器(Filter)
# 以下这种方式导入也可以使用logging，因为logging是包名，导入时会自动执行__init__文件
import logging.handlers
import os
import sys
 
class GetLog(object):
    """单例模式封装日志打印输出到文件和控制台"""
    logger = None
 
    @classmethod
    def get_log(cls, folderName):
        """folderName为日志存放的文件夹名"""
        # logs_path = os.path.dirname(os.path.dirname(__file__)) + '/' + folderName + '/'
        logs_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), folderName + '/')
        if not os.path.exists(logs_path):
            print("ERROR:未找到指定文件夹:%s，请核实！" % folderName)
            sys.exit()
        else:
            if cls.logger is None:
                # 创建日志器，记录器名称通过__name__获取
                cls.logger = logging.getLogger(__name__)
                # 设置日志器级别
                cls.logger.setLevel(logging.INFO)
                # 指定日志文件名为当前执行的py文件名
                filename = str(os.path.basename(sys.argv[0]).split(".")[0]) + '.log'
                # 指定输出的格式和内容
                fmt = '%(asctime)s [%(filename)s-(%(funcName)s:%(lineno)d)] %(levelname)s:%(message)s'
                # 指定日期格式
                datefmt = '%Y-%m-%d %H:%M:%S'
                # 设置格式器
                format = logging.Formatter(fmt, datefmt=datefmt)
 
                # 日志输出到文件，只保存到单一文件，容量会无限大，一般不用，采用时间节点分割(运行时，即使没有日志输出结果，也会创建一个空的日志文件)
                file_handler = logging.FileHandler(logs_path + filename)
                # 日志输出到文件，依据时间切割到不同文件中,when以什么时间单位进行分割（秒、分钟、小时、天、午夜等）
                # interval间隔多久，backupCount保留文件数量，超过数量时间久远的自动被替换掉
                # 例如：以每一晚进行切割日志，保留30个文件
                # file_handler = logging.handlers.TimedRotatingFileHandler(logs_path + filename, when='midnight', interval=1, backupCount=30)
                file_handler.setLevel(logging.INFO)
                file_handler.setFormatter(format)
 
                # 只输出error级别的日志到文件，可以快速定位代码问题(运行时，即使没有日志输出结果，也会创建一个空的日志文件)
                # error_file_handler = logging.handlers.TimedRotatingFileHandler(logs_path + 'error_' + filename, when='midnight', interval=1, backupCount=30)
                # error_file_handler = logging.FileHandler(logs_path + 'error_' + filename)
                # error_file_handler.setLevel(logging.ERROR)
                # error_file_handler.setFormatter(format)
 
                # 使用StreamHandler输出到屏幕
                console = logging.StreamHandler()
                console.setLevel(logging.INFO)
                console.setFormatter(format)
 
                # 添加处理器Handler到日志器中
                cls.logger.addHandler(file_handler)
                # cls.logger.addHandler(error_file_handler)
                cls.logger.addHandler(console)
            return cls.logger
 
if __name__ == "__main__":
    logger = GetLog().get_log('logs')
    logger.info("经营许可证字段值个数:%s", 3)
    logger.info("经营许可证字段值个数:{}".format(333))
    logger.error("错误的日志。。。。")