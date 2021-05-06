
# or
#client = Client(processes=False)  # start local workers as threads
from dask.distributed import Client

if __name__ == '__main__':
    client =Client(process=False)

    def inc(x):
        return x + 1

    def add(x, y):
        return x + y

    a = client.submit(inc, 10)  # calls inc(10) in background thread or process
    b = client.submit(inc, 20)  # calls inc(20) in background thread or process
