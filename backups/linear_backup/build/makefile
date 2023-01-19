CC=gcc
CFLAGS=-O3

SRC_DIR = ../src
BUILD_DIR = .
INSTALL_DIR = ../bin
EXECUTABLE = simple_pm

vpath %.c $(SRC_DIR)
vpath $(EXECUTABLE) $(INSTALL_DIR)

OBJ_NAMES = simple_pm
OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

move : $(EXECUTABLE)
	if [ -e "$(EXECUTABLE)" ]; then mv $(EXECUTABLE) $(INSTALL_DIR);fi 

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $^ -o $@

$(OBJECTS) : %.o : %.c
	$(CC) -o $(EXECUTABLE) $^ $(CFLAGS)

.PHONY : clean clobber
clean:
	-rm -rf *.o
