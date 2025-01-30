import time


def snake_to_camel(input: str) -> str:
    # Can swap out with more sophisticated implementation if needed
    camel_cased = "".join(x.capitalize() for x in input.lower().split("_"))
    if camel_cased:
        return camel_cased[0].lower() + camel_cased[1:]
    else:
        return camel_cased


def time_it(func): 
    '''Decorator that reports the execution time.'''
  
    def wrap(*args, **kwargs): 
        start_time = time.perf_counter()
        result = func(*args, **kwargs) 
        end_time = time.perf_counter()
          
        print(func.__name__, end_time - start_time) 
        return result 
    return wrap 