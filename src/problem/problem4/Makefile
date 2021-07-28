# problem4/Makefile
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
PROG = problem4-aco_acs problem4-aco_as problem4-aco_eas problem4-aco_hcf problem4-aco_mmas problem4-aco_rbas problem4-bd_ig problem4-grasp1 problem4-hsc problem4-urd_ig

.PHONY: all clean debug

all: $(PROG)

clean:
		rm -f *.o *~ core $(PROG)

debug:
		@$(MAKE) OPT+=-D__DEBUG__

problem4-aco_acs:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_acs.o $(LIBS)

problem4-aco_as:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_as.o $(LIBS)

problem4-aco_eas:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_eas.o $(LIBS)

problem4-aco_hcf:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_hcf.o $(LIBS)

problem4-aco_mmas:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_mmas.o $(LIBS)

problem4-aco_rbas:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_rbas.o $(UDIR)/minmax_heap.o $(LIBS)

problem4-bd_ig:		problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/bd_ig.o $(LIBS)

problem4-grasp1:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/grasp1.o $(LIBS)

problem4-hsc:		problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/hsc.o $(LIBS)

problem4-urd_ig:	problem4.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/urd_ig.o $(LIBS)

%.o: %.c
		$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

problem4.o: problem4.c
