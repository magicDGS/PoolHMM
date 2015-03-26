# Copyright (c) 2012 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of LEMON.
#
# LEMON is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# Necessary classes from protable queue from:
# https://github.com/vterron/lemon/blob/9ca6b4b1212228dbd4f69b88aaf88b12952d7d6f/methods.py
# Added by Daniel Gomez-Sanchez

import multiprocessing
import multiprocessing.queues

class SharedCounter(object):
	""" A synchronized shared counter.

	The locking done by multiprocessing.Value ensures that only a single
	process or thread may read or write the in-memory ctypes object. However,
	in order to do n += 1, Python performs a read followed by a write, so a
	second process may read the old value before the new one is written by the
	first process. The solution is to use a multiprocessing.Lock to guarantee
	the atomicity of the modifications to Value.

	This class comes almost entirely from Eli Bendersky's blog:
	http://eli.thegreenplace.net/2012/01/04/shared-counter-with-pythons-multiprocessing/
	"""

	def __init__(self, n = 0):
		self.count = multiprocessing.Value('i', n)

	def increment(self, n = 1):
		""" Increment the counter by n (default = 1) """
		with self.count.get_lock():
			self.count.value += n

	@property
	def value(self):
		""" Return the value of the counter """
		return self.count.value

class Queue(multiprocessing.queues.Queue):
	'''
	A portable implementation of multiprocessing.Queue.

	Because of multithreading / multiprocessing semantics, Queue.qsize() may
	raise the NotImplementedError exception on Unix platforms like Mac OS X
	where sem_getvalue() is not implemented. This subclass addresses this
	problem by using a synchronized shared counter (initialized to zero) and
	increasing / decreasing its value every time the put() and get() methods
	are called, respectively. This not only prevents NotImplementedError from
	being raised, but also allows us to implement a reliable version of both
	qsize() and empty()
	'''

	def __init__(self, *args, **kwargs):
		super(Queue, self).__init__(*args, **kwargs)
		self.size = SharedCounter(0)

	def put(self, *args, **kwargs):
		self.size.increment(1)
		super(Queue, self).put(*args, **kwargs)

	def get(self, *args, **kwargs):
		self.size.increment(-1)
		return super(Queue, self).get(*args, **kwargs)

	def qsize(self):
		'''Reliable implementation of multiprocessing.Queue.qsize()'''
		return self.size.value

	def empty(self):
		'''Reliable implementation of multiprocessing.Queue.empty()'''
		return not self.qsize()