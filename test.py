import os
import datetime
import socket


def get_local_ip():
    hostname = socket.gethostname()
    local_ip = socket.gethostbyname(hostname)
    return local_ip


def system_info(time_now: bool) -> None:
    '''
    Вывод информации о системе
    :param time_now: указиывается, нужно ли вывести текущее время
    :return:
    '''
    print('Windows' if os.name == 'nt' else 'Linux')
    print(os.getlogin())
    if time_now:
        print(datetime.datetime.now().time())

system_info(time_now=True)
print(get_local_ip())
