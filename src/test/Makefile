# test/Makefile
# 
# (C) 2021 Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
# 		   Carlos M. Fonseca <cmfonsec@dei.uc.pt>
# 		   Samuel B. Outeiro <souteiro@student.dei.uc.pt>
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

PROBNAME = default
include macros/default.mk
include macros/$(PROBNAME).mk
include macros/macros.mk

IDIR = ../include
CC = gcc
OPT = -O3 -g
CFLAGS = -I$(IDIR) -Wall -Wextra
LIBS = -lgsl -lgslcblas -lm

.PHONY: all clean debug .FORCE

all: test-$(PROBNAME)

clean:
		rm -f *.o *~ core test-*

debug:
		@$(MAKE) OPT+=-D__DEBUG__

test-$(PROBNAME):	property.o $(PROBPATH)$(PROBOBJ) $(PROBINIT) test.c
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ -I$(PROBPATH) $(LIBS)

property.o:			property.c .FORCE
		$(CC) -c -o $@ $(CFLAGS) $(OPT) $(LIBS) $< $(MACROS)
