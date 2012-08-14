#!/usr/bin/env python

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials
#       provided with the distribution.
# 
#     * Neither the name of Kirk Strauser nor the names of other
#       contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
Very basic parallel processing support

Implements a work-alike of the builtin map() function that distributes
work across many processes.  As it uses Parallel Python to do the
actual multi-processing, code using this must conform to the usual PP
restrictions (arguments must be serializable, etc.)
"""

__author__  = "Kirk Strauser <kirk@strauser.com>"
__version__ = "$Rev: 1139 $"
__date__    = "$Date: 2008-04-16 $"

import time
import __builtin__

try:
    import pp
except:
    raise Exception('Could not import pp')
__STATE = {'server': None}

def ppmap(processes, function, sequence, *sequences):
    """Split the work of 'function' across the given number of
    processes.  Set 'processes' to None to let Parallel Python
    autodetect the number of children to use.

    Although the calling semantics should be identical to
    __builtin__.map (even using __builtin__.map to process
    arguments), it differs in that it returns a generator instead of a
    list.  This enables lazy evaluation of the results so that other
    work can be done while the subprocesses are still running.

    >>> def rangetotal(n): return n, sum(range(n))
    >>> list(map(rangetotal, range(1, 6)))
    [(1, 0), (2, 1), (3, 3), (4, 6), (5, 10)]
    >>> list(ppmap(1, rangetotal, range(1, 6)))
    [(1, 0), (2, 1), (3, 3), (4, 6), (5, 10)]
    """

    # Create a new server if one isn't already initialized
    if not __STATE['server']:
        __STATE['server'] = pp.Server()
    
    def submit(*args):
        """Send a job to the server"""
        return __STATE['server'].submit(function, args, globals=globals())

    # Merge all the passed-in argument lists together.  This is done
    # that way because as with the map() function, at least one list
    # is required but the rest are optional.
    a = [sequence]
    a.extend(sequences)
    available_processes = sum( __STATE['server'].get_active_nodes().values())
    # Set the requested level of multi-processing
    if available_processes < processes:
        __STATE['server'].set_ncpus(processes-available_processes or 'autodetect')

    # First, submit all the jobs.  Then harvest the results as they
    # come available.
    return (subproc() for subproc in __builtin__.map(submit, *a))

if __name__ == '__main__':
    def add(x, y, z):
        """Add three values"""
        return x + y + z

    def busybeaver(x):
        """This can take a while"""
        for num in range(1000000):
            x = x + num
        return x
    # Immediate evaluation example
    start = time.time()
    results = ppmap(None, busybeaver, range(10))
    print 'Time to queue the jobs:', time.time() - start
    start = time.time()
    # Casting the ppmap generator to a list forces each result to be
    # evaluated.  When done immediately after the jobs are submitted,
    # our program twiddles its thumbs while the work is finished.
    print list(results)
    print 'Time to get the results:', time.time() - start

    # Delayed evaluation example
    start = time.time()
    results = ppmap(None, busybeaver, range(10))
    print 'Time to queue the jobs:', time.time() - start
    # In contrast with the above example, this time we're submitting a
    # batch of jobs then going off to do more work while they're
    # processing.  Maybe "time.sleep" isn't the most exciting example,
    # but it illustrates the point that our main program can do work
    # before ppmap() is finished.  Imagine that you're submitting some
    # heavyweight image processing jobs at the beginning of your
    # program, going on to do other stuff like fetching more work to
    # do from a remote server, then coming back later to handle the
    # results.
    time.sleep(5)
    start = time.time()
    print list(results)
    print 'Time to get the first results:', time.time() - start

    # Built-in map example
    print map(add, [1, 2, 3], [4, 5, 6], [7, 8, 9])

    # Trivial ppmap tests
    for i in range(10):
        print '-' * 30
        start = time.time()
        print i, 'adders'
        print ppmap(i, add, [1, 2, 3], [4, 5, 6], [7, 8, 9])
        print 'Iteration time:', time.time() - start

    # Heavier ppmap tests
    for i in range(10):
        print '-' * 30
        start = time.time()
        print i, 'beavers'
        print ppmap(i, busybeaver, range(10))
        print 'Iteration time:', time.time() - start
