# solver/Makefile
#
# (C) 2021 Samuel B. Outeiro <souteiro@student.dei.uc.pt>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

IDIR = ../include
CC = gcc
OPT = -O3 -g
CFLAGS = -I$(IDIR) -Wall -Wextra

.PHONY: all clean debug

all: minmax_heap.o

clean:
		rm -f *.o *~ core

debug:
		@$(MAKE) OPT+=-D__DEBUG__

%.o: %.c
		${CC} -c -o $@ $(CFLAGS) $(OPT) $<

minmax_heap.o: minmax_heap.c
