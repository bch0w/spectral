"""
THIS IS AN UNFINISHED EXAMPLE 8/15 - BC

Examples to show off the different ways to use concurrent.futures.submit() vs.
concurrent.futures.map()

submit: more simple 'scattergun' approach, just submit all jobs without a worry
    about returns, timeouts or exceptions.
map: slightly more complex 'boomerang' approach, submit jobs but expect returns,
    cancel after timeout time, raise exceptions
"""
from concurrent.futures import ProcessPoolExecutor, wait, as_completed


def task(x, break_on=5):
    """Quadratic y=x^2 that breaks on a specific value to show Exceptions"""
    if x == break_on:
        raise Exception
    return x**2


def map_example():
    X = range(0, 11, 1)
    Y = {}
    break_on = 2
    with ProcessPoolExecutor() as executor:
        for x, y in zip(X, executor.map(task, X, break_on, timeout=1)):
            try:
                Y[x] = y
            except Exception:
                Y[x] = "BREAK"
    print(Y)


def submit_example():
    X = range(0, 11, 1)
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(task, x) for x in X]
    done, not_done = wait(futures, timeout=1, return_when="FIRST_EXCEPTION")

    for future in futures:
        print(future.result())  # EXCEPTION

if __name__ == "__main__":
    map_example()
    submit_example()



