# problem2/Makefile
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
PROG = problem2-aco_acs problem2-aco_as problem2-aco_eas problem2-aco_hcf problem2-aco_mmas problem2-aco_rbas problem2-bd_ig problem2-grasp1 problem2-hsc problem2-urd_ig

.PHONY: all clean debug

all: $(PROG)

clean:
		rm -f *.o *~ core $(PROG)

debug:
		@$(MAKE) OPT+=-D__DEBUG__

problem2-aco_acs:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_acs.o $(LIBS)

problem2-aco_as:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_as.o $(LIBS)

problem2-aco_eas:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_eas.o $(LIBS)

problem2-aco_hcf:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_hcf.o $(LIBS)

problem2-aco_mmas:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_mmas.o $(LIBS)

problem2-aco_rbas:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_rbas.o $(UDIR)/minmax_heap.o $(LIBS)

problem2-bd_ig:		problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/bd_ig.o $(LIBS)

problem2-grasp1:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/grasp1.o $(LIBS)

problem2-hsc:		problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/hsc.o $(LIBS)

problem2-urd_ig:	problem2.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/urd_ig.o $(LIBS)

%.o: %.c
		$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

problem2.o: problem2.c
