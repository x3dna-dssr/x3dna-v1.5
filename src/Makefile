### 3DNA Makefile (xiangjun@x3dna.org)

CFLAGS		= -ansi -pedantic -Wall -W -O3
CC		= gcc
LIBS		= -lm

OBJS		= nrutil.o reb_fncs.o ana_fncs.o cmn_fncs.o

HEADERS		= nrutil_3dna.h nrutil.h x3dna.h x3dna_header.h

PROGRAMS	=\
		alc2img \
		analyze \
		anyhelix \
		cehs \
		comb_str \
		ex_str \
		fiber \
		find_pair \
		frame_mol \
		get_part \
		o1p_o2p \
		pdb2img \
		r3d_atom \
		rebuild \
		regular_dna \
		rotate_mol \
		stack2img \
		std_base \
		step_hel

all:	$(PROGRAMS)

# [01]: alc2img
alc2img:	$(OBJS) alc2img.o
		$(CC) $(OBJS) alc2img.o -o alc2img $(LIBS)
		mv alc2img ../bin
alc2img.o:	$(HEADERS)

# [02]: analyze
analyze:	$(OBJS) analyze.o
		$(CC) $(OBJS) analyze.o -o analyze $(LIBS)
		mv analyze ../bin
analyze.o:	$(HEADERS)

# [03]: anyhelix
anyhelix:	$(OBJS) anyhelix.o
		$(CC) $(OBJS) anyhelix.o -o anyhelix $(LIBS)
		mv anyhelix ../bin
anyhelix.o:	$(HEADERS)

# [04]: cehs
cehs:	$(OBJS) cehs.o
		$(CC) $(OBJS) cehs.o -o cehs $(LIBS)
		mv cehs ../bin
cehs.o:	$(HEADERS)

# [05]: comb_str
comb_str:	$(OBJS) comb_str.o
		$(CC) $(OBJS) comb_str.o -o comb_str $(LIBS)
		mv comb_str ../bin
comb_str.o:	$(HEADERS)

# [06]: ex_str
ex_str:	$(OBJS) ex_str.o
		$(CC) $(OBJS) ex_str.o -o ex_str $(LIBS)
		mv ex_str ../bin
ex_str.o:	$(HEADERS)

# [07]: fiber
fiber:	$(OBJS) fiber.o
		$(CC) $(OBJS) fiber.o -o fiber $(LIBS)
		mv fiber ../bin
fiber.o:	$(HEADERS)

# [08]: find_pair
find_pair:	$(OBJS) find_pair.o
		$(CC) $(OBJS) find_pair.o -o find_pair $(LIBS)
		mv find_pair ../bin
find_pair.o:	$(HEADERS)

# [09]: frame_mol
frame_mol:	$(OBJS) frame_mol.o
		$(CC) $(OBJS) frame_mol.o -o frame_mol $(LIBS)
		mv frame_mol ../bin
frame_mol.o:	$(HEADERS)

# [10]: get_part
get_part:	$(OBJS) get_part.o
		$(CC) $(OBJS) get_part.o -o get_part $(LIBS)
		mv get_part ../bin
get_part.o:	$(HEADERS)

# [11]: o1p_o2p
o1p_o2p:	$(OBJS) o1p_o2p.o
		$(CC) $(OBJS) o1p_o2p.o -o o1p_o2p $(LIBS)
		mv o1p_o2p ../bin
o1p_o2p.o:	$(HEADERS)

# [12]: pdb2img
pdb2img:	$(OBJS) pdb2img.o
		$(CC) $(OBJS) pdb2img.o -o pdb2img $(LIBS)
		mv pdb2img ../bin
pdb2img.o:	$(HEADERS)

# [13]: r3d_atom
r3d_atom:	$(OBJS) r3d_atom.o
		$(CC) $(OBJS) r3d_atom.o -o r3d_atom $(LIBS)
		mv r3d_atom ../bin
r3d_atom.o:	$(HEADERS)

# [14]: rebuild
rebuild:	$(OBJS) rebuild.o
		$(CC) $(OBJS) rebuild.o -o rebuild $(LIBS)
		mv rebuild ../bin
rebuild.o:	$(HEADERS)

# [15]: regular_dna
regular_dna:	$(OBJS) regular_dna.o
		$(CC) $(OBJS) regular_dna.o -o regular_dna $(LIBS)
		mv regular_dna ../bin
regular_dna.o:	$(HEADERS)

# [16]: rotate_mol
rotate_mol:	$(OBJS) rotate_mol.o
		$(CC) $(OBJS) rotate_mol.o -o rotate_mol $(LIBS)
		mv rotate_mol ../bin
rotate_mol.o:	$(HEADERS)

# [17]: stack2img
stack2img:	$(OBJS) stack2img.o
		$(CC) $(OBJS) stack2img.o -o stack2img $(LIBS)
		mv stack2img ../bin
stack2img.o:	$(HEADERS)

# [18]: std_base
std_base:	$(OBJS) std_base.o
		$(CC) $(OBJS) std_base.o -o std_base $(LIBS)
		mv std_base ../bin
std_base.o:	$(HEADERS)

# [19]: step_hel
step_hel:	$(OBJS) step_hel.o
		$(CC) $(OBJS) step_hel.o -o step_hel $(LIBS)
		mv step_hel ../bin
step_hel.o:	$(HEADERS)

nrutil.o:	$(HEADERS)
reb_fncs.o:	$(HEADERS)
ana_fncs.o:	$(HEADERS)
cmn_fncs.o:	$(HEADERS)

.PHONY : clean
clean:
	rm *.o
