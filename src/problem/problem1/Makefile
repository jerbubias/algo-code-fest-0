# problem1/Makefile
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

IDIR = ../../include
SDIR = ../../solver
UDIR = ../../util
CC = gcc
OPT = -O3 -g
CFLAGS = -I$(IDIR) -Wall -Wextra
LIBS = -lgsl -lgslcblas -lm
PROG = problem1-aco_acs problem1-aco_as problem1-aco_eas problem1-aco_hcf problem1-aco_mmas problem1-aco_rbas problem1-bd_ig problem1-grasp1 problem1-hsc problem1-urd_ig

.PHONY: all clean debug

all: $(PROG)

clean:
		rm -f *.o *~ core $(PROG)

debug:
		@$(MAKE) OPT+=-D__DEBUG__

problem1-aco_acs:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_acs.o $(LIBS)

problem1-aco_as:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_as.o $(LIBS)

problem1-aco_eas:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_eas.o $(LIBS)

problem1-aco_hcf:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_hcf.o $(LIBS)

problem1-aco_mmas:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_mmas.o $(LIBS)

problem1-aco_rbas:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_rbas.o $(UDIR)/minmax_heap.o $(LIBS)

problem1-bd_ig:		problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/bd_ig.o $(LIBS)

problem1-grasp1:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/grasp1.o $(LIBS)

problem1-hsc:		problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/hsc.o $(LIBS)

problem1-urd_ig:	problem1.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/urd_ig.o $(LIBS)

%.o: %.c
		$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

problem1.o: problem1.c
