import line_profiler

def test():
    print("test")

lp = line_profiler.LineProfiler()
lp.add_function(test)
lp.runctx('test()', locals=locals(), globals=globals())
lp.print_stats()