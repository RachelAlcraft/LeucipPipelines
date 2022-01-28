'''
Taken from https://stackoverflow.com/questions/1316767/how-can-i-explicitly-free-memory-in-python
This is a multithreaded solution to the memory aloocation problem
'''
from dask import delayed  # this module wraps the multithreading
def f(storage, index, chunk_size):  # the processing function
    # read the chunk of size chunk_size starting at index in the file
    # process it using data in storage if needed
    # append data needed for further computations  to storage 
    return storage

partial_result = delayed([])  # put into the delayed() the constructor for your data structure
# I personally use "delayed(nx.Graph())" since I am creating a networkx Graph
chunk_size = 100  # ideally you want this as big as possible while still enabling the computations to fit in memory
for index in range(0, len(file), chunk_size):
    # we indicates to dask that we will want to apply f to the parameters partial_result, index, chunk_size
    partial_result = delayed(f)(partial_result, index, chunk_size)

    # no computations are done yet !
    # dask will spawn a thread to run f(partial_result, index, chunk_size) once we call partial_result.compute()
    # passing the previous "partial_result" variable in the parameters assures a chunk will only be processed after the previous one is done
    # it also allows you to use the results of the processing of the previous chunks in the file if needed

# this launches all the computations
result = partial_result.compute()

# one thread is spawned for each "delayed" one at a time to compute its result
# dask then closes the tread, which solves the memory freeing issue
# the strange performance issue with gc.collect() is also avoided