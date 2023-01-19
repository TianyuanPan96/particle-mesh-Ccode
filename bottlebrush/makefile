CC=gcc
CFLAGS=-O1 -g

ODIR=src/obj

DEPS = src/header.h

INSTALL_DIR = bin
EXECUTABLE = bottlebrush.out

_OBJ = \
		func_math.o \
		func_ran.o \
		func_pm.o \
		func_initialize.o \
		calc_thermo.o \
		mcmoves.o \
		main.o 

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

pm: $(OBJ)
	$(CC) -o $(EXECUTABLE) $^ $(CFLAGS)
	if [ -e "$(EXECUTABLE)" ]; then mv $(EXECUTABLE) $(INSTALL_DIR);fi 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o 
	rm $(INSTALL_DIR)/$(EXECUTABLE)
