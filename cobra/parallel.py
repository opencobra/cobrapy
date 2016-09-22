# Copyright 2013 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

import logging
from multiprocessing import Pool, cpu_count

import six
from six.moves import map
from six.moves import range

from cobra.util import Singleton

logger = logging.getLogger(__name__)


class MultiprocessingView(Singleton):
    """Provides a parallel view (similar to IPython)"""

    def __init__(self, processes=cpu_count(), **kwargs):
        self._processes = processes
        self._kwargs = kwargs
        if not hasattr(self, '_pool'):
            self._pool = None

    @property
    def pool(self):
        if self._pool is None:
            self._pool = Pool(processes=self._processes, **self._kwargs)
        return self._pool

    def __getstate__(self):
        return {
            'processes': self._processes,
            'kwargs': self._kwargs
        }

    def __setstate__(self, state):
        self._kwargs = state['kwargs']
        self._processes = state['processes']

    def map(self, *args, **kwargs):
        return self.pool.map(*args, **kwargs)

    def apply(self, func, *args, **kwargs):
        return self.pool.apply(func, args=args, **kwargs)

    def apply_async(self, func, *args, **kwargs):
        self.pool.apply_async(func, args=args, **kwargs)

    def imap(self, func, *args, **kwargs):
        return self.pool.imap(func, *args, **kwargs)

    def __len__(self):
        return self._processes

    def shutdown(self):
        if self._pool is not None:
            logger.debug('Terminating multiprocessing pool')
            try:
                self._pool.terminate()
            except Exception as e:
                logger.debug('Could not terminate multiprocessing pool.')
                raise e
            finally:
                # TODO: fix the pool
                self._pool = None
        else:
            logger.debug('No multiprocessing pool to shut down.')

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.shutdown()


try:
    import redis
    import six.moves.cPickle as pickle

    class RedisQueue(object):
        """
        Queue with Redis Backend
        http://peter-hoffmann.com/2012/python-simple-queue-redis-queue.html
        """

        MAX_REDIS_LIST_SIZE = 4294967295

        default_host = "localhost"
        default_port = "6379"
        default_db = 0

        def __init__(self, name, maxsize=0, namespace='queue',
                     **connection_args):
            """The default connection parameters are: host='localhost',
            port=6379, db=0

            Parameters
            ----------
            name: str
            maxsize: int
            namespace: str
            connection_args: *kwargs

            """
            if maxsize <= 0:
                maxsize = self.MAX_REDIS_LIST_SIZE
            self._maxsize = maxsize
            self._connection_args = {
                'host': self.default_host,
                'port': self.default_port,
                'db': self.default_db
            }
            for key, val in six.iteritems(connection_args):
                self._connection_args[key] = val
            self._db = redis.Redis(**self._connection_args)
            self._key = '%s:%s' % (namespace, name)

        def __getstate__(self):
            return {
                '_maxsize': self._maxsize,
                '_connection_args': self._connection_args,
                '_key': self._key
            }

        def __setstate__(self, d):
            self._maxsize = d['_maxsize']
            self._connection_args = d['_connection_args']
            self._key = d['_key']
            self._db = redis.Redis(**self._connection_args)

        def __len__(self):
            return self.length

        @property
        def length(self):
            return self._db.llen(self._key)

        def empty(self):
            return self.length() == 0

        def put(self, item):
            """
            Inserts an object in the queue.

            Parameters
            ----------
            item : object
                An object to put in the queue

            """
            if self.length >= self._maxsize:
                raise six.moves.queue.Full

            item = pickle.dumps(item)

            self._db.rpush(self._key, item)

        def put_nowait(self, item):
            """
            Same as put, for backward compatibility with Python Queue.

            See also
            --------
            put

            Parameters
            ----------
            item: object
                An object to put in the queue.

            """
            self.put(item)

        def get(self, block=True, timeout=None):
            """
            Retrieves the next item in the queue.

            Parameters
            ----------
            block: bool, default is True
                If true, the queue is blocked until it an object is
                retrieved reaches the timeout.
            timeout: long
                The timeout (in seconds) that the method should wait for the
                queue to return an item. If block is False, time wil ignored.

            Returns
            -------
            item: object

            """
            if block:
                item = self._db.blpop(self._key, timeout=timeout)
                if item:
                    item = item[1]
                logger.debug("Wait...")
            else:
                item = self._db.lpop(self._key)
                logger.debug("No-wait...")

            logger.debug("Item: %s" % [item])
            if item:
                return pickle.loads(item)
            else:
                raise six.moves.queue.Empty

        def get_nowait(self):
            """Equivalent to get(False)."""
            return self.get(False)

        def __del__(self):
            self._db.delete(self._key)
            self._db.connection_pool.disconnect()

except ImportError:
    pass


class SequentialView(object):
    def map(self, *args, **kwargs):
        return list(map(*args, **kwargs))

    def apply(self, func, *args, **kwargs):
        return func(*args, **kwargs)

    def apply_async(self, func, *args, **kwargs):
        return func(*args, **kwargs)

    def imap(self, func, *args, **kwargs):
        return map(func, *args, **kwargs)

    def __len__(self):
        return 1

    def shutdown(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.shutdown()


if __name__ == '__main__':
    from time import sleep

    class FunctionObject(object):
        """docstring for FunctionObject"""

        def __init__(self):
            super(FunctionObject, self).__init__()

        def __call__(self, arg):
            sleep(.1)
            return arg ** 2

    view = MultiprocessingView(processes=4)
    print(view.map(FunctionObject(), list(range(0, 100))))
