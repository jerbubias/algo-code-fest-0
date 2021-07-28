# problem3/Makefile
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
PROG = problem3-aco_acs problem3-aco_as problem3-aco_eas problem3-aco_hcf problem3-aco_mmas problem3-aco_rbas problem3-bd_ig problem3-grasp1 problem3-hsc problem3-urd_ig

.PHONY: all clean debug

all: $(PROG)

clean:
		rm -f *.o *~ core $(PROG)

debug:
		@$(MAKE) OPT+=-D__DEBUG__

problem3-aco_acs:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_acs.o $(LIBS)

problem3-aco_as:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_as.o $(LIBS)

problem3-aco_eas:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_eas.o $(LIBS)

problem3-aco_hcf:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_hcf.o $(LIBS)

problem3-aco_mmas:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_mmas.o $(LIBS)

problem3-aco_rbas:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_rbas.o $(UDIR)/minmax_heap.o $(LIBS)

problem3-bd_ig:		problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/bd_ig.o $(LIBS)

problem3-grasp1:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/grasp1.o $(LIBS)

problem3-hsc:		problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/hsc.o $(LIBS)

problem3-urd_ig:	problem3.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/urd_ig.o $(LIBS)

%.o: %.c
		$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

problem3.o: problem3.c
