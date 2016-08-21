
Some of my VISTA Hemisphere Python code. Use at your own risk. I have been 
migrating legacy IDL code to IDL and will put the Python versions here. My
Python style is evolving from FORTRAN/IDL to become Python 3.x and PEP 08 
compliant.

My code development style is usually as follows:

1.  write monolithic code to deal with short term research problem

2.  when it become unwieldy I refactor into functions and use ```__main__``` to run 
the original i.e.

def task1():

```
if __name__ == '__main__':

     your code here
```



You will probably need librgm to run this code. 

To run with dependencies:

1.  do what you usually do depending on your own preferences OR

2.  dump all you need into the same directory OR

3.  use sys.path.append for functions not in PYTHONPATH
     e.g. 

```sys.path.append('/home/rgm/soft/python/lib/')```

Some of the local data dependencies are stored in the config file vhs.ini. 
