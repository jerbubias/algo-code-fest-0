# mkp/Makefile
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
PROG = mkp-aco_acs mkp-aco_as mkp-aco_eas mkp-aco_hcf mkp-aco_mmas mkp-aco_rbas mkp-bd_ig mkp-d_df_bnb mkp-grasp1 mkp-grasp2 mkp-hsc mkp-p_df_bnb mkp-urd_ig

.PHONY: all clean debug

all: $(PROG)

clean:
		rm -f *.o *~ core $(PROG)

debug:
		@$(MAKE) OPT+=-D__DEBUG__

mkp-aco_acs:	mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_acs.o $(LIBS)

mkp-aco_as:		mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_as.o $(LIBS)

mkp-aco_eas:	mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_eas.o $(LIBS)

mkp-aco_hcf:	mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_hcf.o $(LIBS)

mkp-aco_mmas:	mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_mmas.o $(LIBS)

mkp-aco_rbas:	mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/aco_rbas.o $(UDIR)/minmax_heap.o $(LIBS)

mkp-bd_ig:		mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/bd_ig.o $(LIBS)

mkp-d_df_bnb:	mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/d_df_bnb.o $(LIBS)

mkp-grasp1:		mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/grasp1.o $(LIBS)

mkp-grasp2:		mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/grasp2.o $(LIBS)

mkp-hsc:		mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/hsc.o $(LIBS)

mkp-p_df_bnb:	mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/p_df_bnb.o $(LIBS)

mkp-urd_ig:		mkp.o
		$(CC) -o $@ $(CFLAGS) $(OPT) $^ $(SDIR)/urd_ig.o $(LIBS)

%.o: %.c
		$(CC) -c -o $@ $(CFLAGS) $(OPT) $<

mkp.o: mkp.c
